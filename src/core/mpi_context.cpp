#include "core/mpi_context.hpp"

#include <array>
#include <stdexcept>

MpiContext::MpiContext(MPI_Comm comm) : comm_(comm) {
    check(MPI_Comm_set_errhandler(comm_, MPI_ERRORS_RETURN), "MPI_Comm_set_errhandler");
    check(MPI_Comm_rank(comm_, &rank_), "MPI_Comm_rank");
    check(MPI_Comm_size(comm_, &size_), "MPI_Comm_size");
}

void MpiContext::check(int error, const std::string& operation) const {
    if (error == MPI_SUCCESS) return;
    std::array<char, MPI_MAX_ERROR_STRING> message{};
    int length = 0;
    MPI_Error_string(error, message.data(), &length);
    throw std::runtime_error(operation + " failed on rank " + std::to_string(rank_) + ": " +
                             std::string(message.data(), static_cast<std::size_t>(length)));
}

void MpiContext::sendrecv_doubles(const double* send, int send_count, int peer, int send_tag,
                                  double* recv, int recv_count, int recv_tag,
                                  const char* operation) const {
    const int error = MPI_Sendrecv(send, send_count, MPI_DOUBLE, peer, send_tag,
                                   recv, recv_count, MPI_DOUBLE, peer, recv_tag,
                                   comm_, MPI_STATUS_IGNORE);
    check(error, std::string(operation) + " peer=" + std::to_string(peer) +
                     " send_tag=" + std::to_string(send_tag) +
                     " recv_tag=" + std::to_string(recv_tag));
}

void MpiContext::alltoallv_doubles(const double* send,
                                   const std::vector<int>& send_counts,
                                   const std::vector<int>& send_displs,
                                   double* recv,
                                   const std::vector<int>& recv_counts,
                                   const std::vector<int>& recv_displs,
                                   const char* operation) const {
    const auto communicator_size = static_cast<std::size_t>(size_);
    if (send_counts.size() != communicator_size || send_displs.size() != communicator_size ||
        recv_counts.size() != communicator_size || recv_displs.size() != communicator_size) {
        throw std::invalid_argument(std::string(operation) +
                                    " count/displacement size does not match communicator");
    }
    check(MPI_Alltoallv(send, send_counts.data(), send_displs.data(), MPI_DOUBLE,
                        recv, recv_counts.data(), recv_displs.data(), MPI_DOUBLE, comm_),
          operation);
}

void MpiContext::barrier(const char* operation) const {
    check(MPI_Barrier(comm_), operation);
}

double MpiContext::reduce_sum(double local, int root, const char* operation) const {
    double global = 0.0;
    check(MPI_Reduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, root, comm_), operation);
    return global;
}

void MpiContext::broadcast_doubles(double* values, int count, int root,
                                   const char* operation) const {
    check(MPI_Bcast(values, count, MPI_DOUBLE, root, comm_), operation);
}

int MpiContext::allreduce_min(int local, const char* operation) const {
    int global = 0;
    check(MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MIN, comm_), operation);
    return global;
}

std::vector<std::uint64_t> MpiContext::allgather_u64(std::uint64_t local,
                                                     const char* operation) const {
    std::vector<std::uint64_t> values(static_cast<std::size_t>(size_));
    check(MPI_Allgather(&local, 1, MPI_UINT64_T, values.data(), 1, MPI_UINT64_T, comm_), operation);
    return values;
}

std::vector<std::int64_t> MpiContext::allgather_i64(std::int64_t local,
                                                    const char* operation) const {
    std::vector<std::int64_t> values(static_cast<std::size_t>(size_));
    check(MPI_Allgather(&local, 1, MPI_INT64_T, values.data(), 1, MPI_INT64_T, comm_), operation);
    return values;
}
