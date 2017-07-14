import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, lfilter, cheby2, medfilt

def moving_average(a, n=3):
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def bandpass_coefficients(lowcut, highcut, fs, order=10):
	nyq = 0.5 * fs
	low = lowcut / nyq
	high = highcut / nyq
	b, a = cheby2(order, 50, [low,high], btype='bandpass')
	return b, a

def bandpass(data, low, high, rate, order=None):
	bpf = bandpass_coefficients(low, high, rate, order)
	return lfilter(bpf[0], bpf[1], data)

def band_powers(data, fs, freqs, pulse_duration,
		sample_duration, bandwidth=100, debug=False):
	pulse_length = int(fs * pulse_duration)
	sample_length = int(fs * sample_duration)
	tones = []
	powers = {}
	for freq in freqs:
		avg_powers = []
		for channel in xrange(data.shape[1]):
			band = bandpass(data[:,channel], freq - bandwidth, freq + bandwidth, fs, 3)
			tone_start = np.argmax(moving_average(band, pulse_length)) - pulse_length
			if channel == 0:
				tones.append(tone_start)
			avg_powers.append(
				(1.0 / pulse_length) * np.sum(
					band[tone_start:tone_start+sample_length] ** 2
				)
			)
		powers[freq] = np.mean(avg_powers)
	if debug:
		plt.plot(np.array(tones) / fs, freqs, 'bx')
		plt.specgram(data[:,0], Fs=fs, NFFT=1024)
		plt.show()
	return powers