#pragma once

#include "core/params.hpp"

inline Params stage3_params() {
    return Params::load(CHANNEL_DNS_SOURCE_DIR "/examples");
}
