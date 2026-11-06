#include "core/cpu_device.hpp"

#include <condition_variable>
#include <deque>
#include <exception>
#include <mutex>
#include <thread>

// Private thread-pool backend. No scheduler type or factory is exposed to callers.
class CpuDevice::Impl {
  public:
    explicit Impl(int threads);
    ~Impl();
    Impl(const Impl &) = delete;
    Impl &operator=(const Impl &) = delete;
    void enqueue(std::function<void()> task);
    void sync();
    void require_owner() const;
    void enqueue_batch(std::unique_ptr<Batch> batch);

  private:
    struct SingleTask;
    struct Task {
        Batch *batch;
        std::size_t lane;
    };
    Task take_task(); // mutex_ must be held.
    void worker();
    void execute(Task task);
    void stop();
    const std::thread::id owner_;
    std::vector<std::thread> workers_;
    std::mutex mutex_;
    std::condition_variable ready_, done_;
    std::deque<std::unique_ptr<Batch>> queue_;
    std::exception_ptr error_;
    std::size_t pending_ = 0;
    bool finishing_batch_ = false;
    bool stopping_ = false;
};

thread_local bool CpuDevice::executing_ = false;

struct CpuDevice::Impl::SingleTask final : Batch {
    explicit SingleTask(std::function<void()> function) : Batch(1), function(std::move(function)) {}
    void run(std::size_t) override { function(); }
    std::function<void()> function;
};

CpuDevice::Impl::Impl(int threads) : owner_(std::this_thread::get_id()) {
    if (threads < 1)
        throw std::invalid_argument("thread count must be positive");
    workers_.reserve(static_cast<std::size_t>(threads));
    try {
        for (int i = 0; i < threads; ++i)
            workers_.emplace_back([this] { worker(); });
    } catch (...) {
        stop();
        throw;
    }
}

CpuDevice::Impl::~Impl() {
    try {
        sync();
    } catch (...) {
    } // Explicit sync() is required to observe task errors.
    stop();
}

void CpuDevice::Impl::require_owner() const {
    if (executing_ || std::this_thread::get_id() != owner_)
        throw std::logic_error(
            "queue submission and sync must originate on the owner outside a task");
}

void CpuDevice::Impl::stop() {
    {
        std::lock_guard lock(mutex_);
        stopping_ = true;
    }
    ready_.notify_all();
    for (auto &worker : workers_)
        worker.join();
}

void CpuDevice::Impl::enqueue(std::function<void()> task) {
    require_owner();
    if (!task)
        throw std::invalid_argument("cannot enqueue an empty task");
    enqueue_batch(std::make_unique<SingleTask>(std::move(task)));
}

void CpuDevice::Impl::enqueue_batch(std::unique_ptr<Batch> batch) {
    const auto chunks = batch->chunks;
    bool activate;
    {
        std::lock_guard lock(mutex_);
        // Publish only after all callable copies have been constructed successfully.
        activate = queue_.empty() && !finishing_batch_;
        queue_.push_back(std::move(batch));
        pending_ += chunks;
    }
    if (!activate)
        return; // A preceding batch will wake workers when this one becomes executable.
    if (chunks == 1)
        ready_.notify_one();
    else
        ready_.notify_all();
}

CpuDevice::Impl::Task CpuDevice::Impl::take_task() {
    auto *batch = queue_.front().get();
    const auto lane = batch->next++;
    return {batch, lane};
}

void CpuDevice::Impl::execute(Task task) {
    executing_ = true;
    try {
        task.batch->run(task.lane);
    } catch (...) {
        std::lock_guard lock(mutex_);
        if (!error_)
            error_ = std::current_exception();
    }
    std::unique_lock lock(mutex_);
    if (--task.batch->remaining != 0) {
        executing_ = false;
        --pending_;
        return;
    }
    auto completed = std::move(queue_.front());
    queue_.pop_front();
    finishing_batch_ = true;
    lock.unlock();
    // Destruct user captures outside the mutex, before releasing the next batch.
    completed.reset();
    executing_ = false;
    lock.lock();
    finishing_batch_ = false;
    if (--pending_ == 0)
        done_.notify_one();
    // This worker can run a following single-chunk command directly. Wake peers
    // only when the next batch exposes parallel work, or when shutting down.
    const bool wake_workers = stopping_ || (!queue_.empty() && queue_.front()->chunks > 1);
    lock.unlock();
    if (wake_workers)
        ready_.notify_all();
}

void CpuDevice::Impl::worker() {
    std::unique_lock lock(mutex_);
    for (;;) {
        ready_.wait(lock, [&] {
            return (stopping_ && pending_ == 0) ||
                   (!finishing_batch_ && !queue_.empty() &&
                    queue_.front()->next < queue_.front()->chunks);
        });
        if (stopping_ && pending_ == 0)
            return;
        auto task = take_task();
        lock.unlock();
        execute(std::move(task));
        lock.lock();
    }
}

void CpuDevice::Impl::sync() {
    require_owner();
    std::unique_lock lock(mutex_);
    done_.wait(lock, [&] { return pending_ == 0; });
    auto error = std::exchange(error_, nullptr);
    lock.unlock();
    if (error)
        std::rethrow_exception(error);
}

CpuDevice::CpuDevice(int threads)
    : threads_(threads > 0 ? threads : static_cast<int>(std::max(1u, std::thread::hardware_concurrency()))),
      impl_(std::make_unique<Impl>(threads_)) {}
CpuDevice::~CpuDevice() = default;
void CpuDevice::require_owner() const { impl_->require_owner(); }
void CpuDevice::enqueue(std::function<void()> task) { impl_->enqueue(std::move(task)); }
void CpuDevice::enqueue_batch(std::unique_ptr<Batch> batch) { impl_->enqueue_batch(std::move(batch)); }
void CpuDevice::sync() { impl_->sync(); }

