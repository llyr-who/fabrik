#include "Timer.h"
#include <time.h>

#define OONPS 1E-9;  // ONE_OVER_NANOSECONDS_PER_SECOND

GameTimer::GameTimer() {
    deltaTime = -1.0;
    clock_gettime(CLOCK_REALTIME, &baseTime);
    pausedTime = {0, 0};
    prevTime = {0, 0};
    currTime = {0, 0};
    stopped = false;
}

float GameTimer::GameTime() {
    if (stopped) {
        return (stopTime.tv_sec - pausedTime.tv_sec) - baseTime.tv_sec +
               ((stopTime.tv_nsec - pausedTime.tv_nsec) - baseTime.tv_nsec) *
                   OONPS;
    } else {
        return (currTime.tv_sec - baseTime.tv_sec) +
               (currTime.tv_nsec - baseTime.tv_nsec) * OONPS;
    }
}

void GameTimer::Reset() {
    deltaTime = -1.0;
    clock_gettime(CLOCK_REALTIME, &baseTime);
    pausedTime = {0, 0};
    prevTime = {0, 0};
    currTime = {0, 0};
    stopped = false;
}

void GameTimer::Tick() {
    if (stopped) {
        deltaTime = 0.0;
        return;
    }
    clock_gettime(CLOCK_REALTIME, &currTime);
    deltaTime = (currTime.tv_sec - prevTime.tv_sec) +
                (currTime.tv_nsec - prevTime.tv_nsec) * OONPS;
    prevTime = currTime;
}

float GameTimer::DeltaTime() { return (float)deltaTime; }

void GameTimer::Start() {
    struct timespec start;
    clock_gettime(CLOCK_REALTIME, &start);

    // Accumulate the time elapsed between stop and start pairs.
    //
    //                     |<-------d------->|
    // ----*---------------*-----------------*------------> time
    //  mBaseTime       mStopTime        startTime

    if (stopped) {
        pausedTime.tv_sec += (start.tv_sec - stopTime.tv_sec);
        pausedTime.tv_nsec += (start.tv_nsec - stopTime.tv_nsec);

        prevTime = start;
        stopTime.tv_sec = 0;
        stopTime.tv_nsec = 0;
        stopped = false;
    }
}

void GameTimer::Stop() {
    if (!stopped) {
        struct timespec current;
        clock_gettime(CLOCK_REALTIME, &current);

        stopTime = currTime;
        stopped = true;
    }
}

