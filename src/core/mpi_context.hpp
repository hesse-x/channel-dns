#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include <mpi.h>

class MpiContext {
public:
    explicit MpiContext(MPI_Comm comm);

    MPI_Comm comm() const { return comm_; }
    int rank() const { return rank_; }
    int size() const { return size_; }

    void sendrecv_doubles(const double* send, int send_count, int peer, int send_tag,
                          double* recv, int recv_count, int recv_tag,
                          const char* operation) const;
    void alltoallv_doubles(const double* send,
                           const std::vector<int>& send_counts,
                           const std::vector<int>& send_displs,
                           double* recv,
                           const std::vector<int>& recv_counts,
                           const std::vector<int>& recv_displs,
                           const char* operation) const;
    void barrier(const char* operation) const;
    double reduce_sum(double local, int root, const char* operation) const;
    void broadcast_doubles(double* values, int count, int root, const char* operation) const;
    int allreduce_min(int local, const char* operation) const;
    std::vector<std::uint64_t> allgather_u64(std::uint64_t local,
                                             const char* operation) const;
    std::vector<std::int64_t> allgather_i64(std::int64_t local,
                                            const char* operation) const;

private:
    void check(int error, const std::string& operation) const;

    MPI_Comm comm_;
    int rank_ = 0;
    int size_ = 1;
};
