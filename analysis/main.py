import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import spline
from scikits.audiolab import wavread

from segment import band_powers

SAMPLE_DURATIONS = [d * 0.001 for d in xrange(10, 65)]
FREQS = [25000, 30000, 35000, 40000, 45000]
PULSE_DURATION = 0.06
SWEEP_DURATION = 0.2

def normalize_angle(angle):
	while not -180 <= angle < 180:
		angle += 360 * (1 if angle < 0 else -1)
	return angle * (np.pi / 180)

# References is a dict mapping experiment csv => folder
def get_files(references):
	readings = []
	azimuth_files = {}
	for ref in references:
		path_template = references[ref]
		with open(ref,'r') as csvfile:
			reader = csv.reader(csvfile)
			reader.next()
			for row in reader:
				azimuth = int(row[0])
				files = [path_template % (e.zfill(4)) for e in row[1:]]
				if not azimuth_files.get(azimuth):
					azimuth_files[azimuth] = []
				azimuth_files[azimuth] += files
	return azimuth_files

def extract_powers(azimuth_files, sample_duration):
	azimuth_by_freq = {}
	for f in FREQS:
		azimuth_by_freq[f] = {}
	for a in sorted(azimuth_files):
		for file in azimuth_files[a]:
			data, fs, enc = wavread(file)
			print "Calculating azimuth powers", a, file
			powers = band_powers(data, fs, FREQS,
				PULSE_DURATION, sample_duration, SWEEP_DURATION, debug=False)
			for f in powers:
				if not azimuth_by_freq[f].get(a):
					azimuth_by_freq[f][a] = []
				azimuth_by_freq[f][a].append(powers[f])
	return azimuth_by_freq

def normalize_band(azimuths):
	assert azimuths[0]
	for a in azimuths:
		azimuths[a] = np.median(azimuths[a])
	db_azimuths = {}
	for a in azimuths:
		db_azimuths[a] = 10 * np.log10(azimuths[a] / azimuths[0])
	return db_azimuths

def normalize(azimuth_powers):
	db_azimuths_freqs = {}
	for freq in azimuth_powers:
		db_azimuths_freqs[freq] = normalize_band(azimuth_powers[freq])
	return db_azimuths_freqs

def plot_gain_pattern(azimuth_powers, db_azimuths_freqs, filename):
	ax = plt.subplot(111, projection='polar')
	for freq in sorted(db_azimuths_freqs):
		rel = 10 * np.log10(azimuth_powers[freq][0] / azimuth_powers[FREQS[0]][0])
		powers_dict = db_azimuths_freqs[freq]
		azimuths = np.array(map(normalize_angle,sorted(powers_dict.keys())))
		# Aaaaaaah --- we need azimuths to be monotonic for splice to work
		slice_index = - (np.argmax(np.diff(azimuths) < 0) + 1)
		azimuths = np.roll(azimuths, slice_index)
		powers = np.roll(
			np.array(
				[powers_dict[p] + rel for p in sorted(powers_dict.keys())]
			),
			slice_index
		)
		print azimuths, powers
		azimuths_smooth = np.linspace(azimuths.min(), azimuths.max(), 1e3)
		powers_smooth = spline(azimuths, powers, azimuths_smooth)
		ax.plot(azimuths_smooth, powers_smooth, '-', label=str(freq))
		ax.set_theta_zero_location('N')
		ax.set_ylim(-70, 40)
		ax.set_xticks(np.arange(0, 360, 10) * (np.pi / 180))
		ax.set_yticks(np.arange(-70, 40, 10))
		# ax.set_rlabel_position(-np.pi / 2)
		# ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
		#           ncol=len(FREQS), mode="expand", borderaxespad=0, prop={'size': 6})
	plt.savefig(filename)
	plt.close()

def analysis(references, sample_duration):
	azimuth_files = get_files(references)
	azimuth_powers = extract_powers(azimuth_files, sample_duration)
	print azimuth_powers
	db_azimuths_freqs = normalize(azimuth_powers)
	plot_gain_pattern(
		azimuth_powers,
		db_azimuths_freqs, '../output/okay-%s.png' % (sample_duration)
	)

if __name__ == "__main__":
	references = {
		"../samples/2017-7-12/azimuth-2.csv": "../samples/2017-7-12/azimuth-2/DR0000_%s.wav",
		"../samples/2017-7-12/azimuth-1-mini.csv": "../samples/2017-7-12/azimuth-1/DR0000_%s.wav",
	}
	for d in SAMPLE_DURATIONS:
		analysis(references, d)