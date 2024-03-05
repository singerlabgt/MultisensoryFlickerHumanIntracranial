#adapted from Jacques Lerousseau
#this script expects the following input variables with those names:
# sfreq: the sampling frequency
# trial_duration
# prestim_duration:
# poststim_duration:
# condition_freq:
# EVOKED:


#import required packages
import sys
import os
import numpy as np
import mne
import scipy.signal
import scipy.io.wavfile
import re

#set times variable:
times = np.arange(-prestim_duration, trial_duration+poststim_duration, 1./sfreq) #timepoints of evoked potential


n_channels, n_times = EVOKED.shape

# Filter EVOKED
EVOKED_FILTERED = mne.filter.filter_data(EVOKED, sfreq, condition_freq-0.5, condition_freq+0.5) #filtered signal +/-1Hz from condition frequency

EVOKED_HIL = np.abs(scipy.signal.hilbert(EVOKED_FILTERED, axis=-1)) #take Hilbert of filtered signal

# Z-score compare to pre-stimulus activity
pre_stim = EVOKED_HIL[:, times <= 0.]
mu = np.mean(pre_stim, axis=-1, keepdims=True)
sd = np.std(pre_stim, axis=-1, keepdims=True)

# Z transform
EVOKED_HIL_Z = (EVOKED_HIL - mu) / sd

# Define threshold; chosen such that no channel/vertice shows activity before stimulus onset
#thresh = np.max(EVOKED_HIL_Z[:, times <= 0.])
thresh=3 #setting to 3

# Find onsets, i.e. the time at which stimulus response is considered to start
ONSETS = np.full((n_channels), np.nan)
for i in np.arange(n_channels):
    _onset = np.where(EVOKED_HIL_Z[i] > thresh)[0]
    if _onset.size == 0: continue
    ONSETS[i] = times[_onset[0]]

# Remove onsets < 0
ONSETS[ONSETS <= 0.] = np.nan

times_rec = np.arange(-15, 15, 1./sfreq)
EVOKED_HIL_REC = np.full((n_channels, times_rec.size), np.nan)
for i, onset in enumerate(ONSETS):
    if np.isnan(onset): continue
    EVOKED_HIL_REC[i, np.logical_and(times_rec >= -1 - onset, times_rec < 11 - onset)] = EVOKED_HIL_Z[i]

#bin_onsets = np.arange(-condition_freq * (1/condition_freq), 11., 1/condition_freq) #EQUIVALENT TO TIMEPOINTS OF EACH ONSET?
#BELOW IS CORRECT? SO DEALS WITH 5.5HZ CONDITIONS
bin_onsets=np.arange(-np.ceil(condition_freq) * (1/condition_freq), 11., 1/condition_freq)
freq = condition_freq #FREQUENCY OF STIM?

EVOKED_HIL_REC_BIN = np.full((bin_onsets.size, n_channels), np.nan)

for j, bin_onset in enumerate(bin_onsets):

    # Find index of time
    tmin = bin_onset
    tmax = bin_onset + 1./freq
    t_idx = np.logical_and(times_rec >= tmin, times_rec < tmax)

    # Average in this time bins
    EVOKED_HIL_REC_BIN[j] = np.mean(EVOKED_HIL_REC[:, t_idx], axis=-1) #averages the z-scored hilbert for each time bin

# Define threshold
threshold = np.nanmax(EVOKED_HIL_REC_BIN[:int(condition_freq)]) #chosen such that no channel/vertice shows activity before response onset CORRECT?
threshold

GOOD_CHANNELS = np.ones((n_channels), dtype=np.bool)

# Check if channel has an onset value.
GOOD_CHANNELS[np.isnan(ONSETS)] = False

# Check if channel is activated on the first bin.
for j in range(n_channels):
    if EVOKED_HIL_REC_BIN[int(condition_freq), j] < threshold: #CORRECT?
        GOOD_CHANNELS[j] = False

# Length of enveloppe
def countConsecutive(X):
    return np.diff(np.where(np.concatenate(([X[0]], X[:-1] != X[1:], [True])))[0])[::2]

N_CYCLES = np.full((n_channels), np.nan)
T = EVOKED_HIL_REC_BIN > threshold
for j in np.arange(n_channels):
    if countConsecutive(T[:, j]).size > 0:
        N_CYCLES[j] = np.max(countConsecutive(T[:, j]))


# Build cumulative distribution
CUM_DISTR = np.zeros((bin_onsets.size))

for i in range(n_channels):
    if np.isnan(N_CYCLES[i]): continue
    CUM_DISTR[int(N_CYCLES[i])] += 1

# Cumulative
CUM_DISTR = np.cumsum(CUM_DISTR)


#variables that should be returned:
#EVOKED_FILTERED
#EVOKED_HIL
#EVOKED_HIL_Z
#thresh
#EVOKED_HIL_REC
#EVOKED_HIL_REC_BIN
#threshold
#GOOD_CHANNELS
#N_CYCLES
#ONSETS
