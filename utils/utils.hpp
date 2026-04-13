#pragma once

// Utility functions for time conversion and other time related operations
namespace TimeUtils {
    [[maybe_unused]] static int convertRawTimeToPhysical(int raw_time, int bcid) {
        return raw_time != 0 ? bcid * 30 + raw_time - 1 : -1;
    }
}