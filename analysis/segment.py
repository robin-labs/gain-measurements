import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft
from scipy.signal import butter, lfilter, freqz
import peakutils

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float) ** 2
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def butter_highpass(cutoff, fs, order=3):
    nyq = 0.5 * fs
    b, a = butter(order, cutoff / nyq, btype='high')
    return b, a

def highpass(data, cutoff, fs, order=3):
    b, a = butter_highpass(cutoff, fs, order)
    return lfilter(b, a, data)

def butter_bandpass(lowcut, highcut, fs, order=3):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def bandpass(data, lowcut, highcut, fs, order=3):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    return lfilter(b, a, data)

def band_powers_many(data, n_iterations, fs, freqs, pulse_duration,
        sample_duration, sweep_duration, marker_distance_s=3, find_bandwidth=500, total_bandwidth=5000, debug=False):
	pulse_length = int(fs * pulse_duration)
	MARKER_FREQ = 1e4
	marker_band = bandpass(data[:,0], MARKER_FREQ - find_bandwidth / 2, MARKER_FREQ + find_bandwidth / 2, fs)
	marker_powers = marker_band ** 2
	marker_peaks = peakutils.indexes(marker_powers, 0.5, fs * marker_distance_s)
	if len(marker_peaks) != n_iterations:
		raise Error("Unexpected number of iterations detected: % instead of %s" % (len(marker_peaks, n_iterations)))
	sections = np.split(data, marker_peaks)[:-1]
	powers = []
	for section in sections:
		print "Computing section..."
		# For each section, find the 50ms 15k pulse and cut everything off before it
		# It's clumsy, but it works
		LOWEST_FREQ = 1.5e4
		find_band = bandpass(section[:,0], LOWEST_FREQ - 100, LOWEST_FREQ + 100, fs)
		start_lowest = np.argmax(moving_average(find_band, pulse_length))
		section = section[0:start_lowest]
		powers.append(band_powers(section, fs, freqs, pulse_duration, sample_duration, 
				sweep_duration, find_bandwidth, total_bandwidth, debug))
	avg_powers_across_samples = {}
	for p in sorted(powers[0]):
		collected = [sample_power[p] for sample_power in powers]
		print "f = %s" % p, "mean = %s" % np.mean(collected), \
			"median = %s" % np.median(collected), "std = %s" % np.std(collected)
		avg_powers_across_samples[p] = np.mean(collected)
	return avg_powers_across_samples


def band_powers(data, fs, freqs, pulse_duration,
        sample_duration, sweep_duration, find_bandwidth=400, total_bandwidth=5000, debug=False):
    channels = xrange(data.shape[1])
    pulse_length = int(fs * pulse_duration)
    sample_length = int(fs * sample_duration)
    tones = []
    powers = {}
    for channel in channels:
        data[:,channel] = highpass(data[:,channel], 1e4, fs, 5)
    for freq in freqs:
        avg_powers = []
        for channel in channels:
            band = bandpass(data[:,channel], freq - find_bandwidth / 2, freq + find_bandwidth / 2, fs)
            tone_start = np.argmax(moving_average(band, pulse_length))
            if channel == 0:
                tones.append(tone_start)
            band = bandpass(data[:,channel], freq - total_bandwidth / 2, freq + total_bandwidth / 2, fs)
            power = np.sum(band[tone_start:tone_start+sample_length] ** 2)
            avg_powers.append(power)
            # data = data[tone_start:]
        powers[freq] = np.mean(avg_powers)
    if debug:
        plt.plot(np.array(tones, dtype=np.float) / fs, freqs, 'bx')
        plt.specgram(data[:,0], Fs=fs, NFFT=1024)
        plt.show()
    return powers