package org.usadellab.trimmomatic.trim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 * Long-read chimera splitting using k-mer-seeded Hamming-distance matching.
 *
 * <p>Scans the interior of each read (excluding the terminal zones at both ends,
 * which are handled by {@link LongReadClipTrimmer}) for adapter sequences. When
 * one or more internal adapter hits are found the read is split into fragments;
 * the adapter bases are discarded. Each fragment is passed independently through
 * all subsequent trimming steps in the pipeline.
 *
 * <p><b>Single-end mode only.</b> Invoking this step in paired-end mode raises
 * a {@link RuntimeException} immediately, because no current long-read platform
 * (ONT, PacBio) produces paired-end data.
 *
 * <p>Recommended pipeline order:
 * <pre>
 *   LONGREADSPLIT:adapters.fa:0.10  LONGREADCLIP:adapters.fa:0.10  MINLEN:200
 * </pre>
 *
 * <p>Step syntax:
 * <pre>
 *   LONGREADSPLIT:&lt;fasta&gt;:&lt;maxErrorRate&gt;[:&lt;minOverlap&gt;[:&lt;minFragmentLength&gt;]]
 * </pre>
 *
 * <p>Parameters:
 * <ul>
 *   <li>{@code fasta} – FASTA file of adapter sequences. Both the forward
 *       orientation and its reverse complement are loaded automatically.</li>
 *   <li>{@code maxErrorRate} – maximum fraction of mismatches allowed in a
 *       confirmed adapter match (e.g. {@code 0.10} = 1 mismatch per 10 bp).
 *       N bases in either the read or the adapter act as wildcards and are
 *       never counted as mismatches.</li>
 *   <li>{@code minOverlap} – minimum number of overlapping bases required to
 *       call an adapter match [default 10]. Also defines the <em>terminal
 *       zone</em>: matches whose start position falls within this many bases
 *       of either end are ignored (handled by LONGREADCLIP).</li>
 *   <li>{@code minFragmentLength} – fragments shorter than this after splitting
 *       are discarded [default 100].</li>
 * </ul>
 *
 * <p>Split fragments are renamed by appending {@code /splitNofM} to the
 * original read name, where N is the 1-based fragment index and M is the
 * total number of fragments produced from that read.
 */
public class LongReadSplitTrimmer implements Trimmer {

	private static final int DEFAULT_MIN_OVERLAP         = 10;
	private static final int DEFAULT_MIN_FRAGMENT_LENGTH = 100;
	private static final int KMER_SIZE                   = 8;

	private final List<String>              adapters;
	/** k-mer → list of {adapterIndex, positionInAdapter} */
	private final Map<String, List<int[]>>  kmerIndex;
	private final float maxErrorRate;
	private final int   minOverlap;
	private final int   minFragmentLength;

	// ------------------------------------------------------------------
	// Construction

	public LongReadSplitTrimmer(String args) throws IOException {
		String[] parts = args.split(":");
		if (parts.length < 2) {
			throw new IllegalArgumentException(
					"LONGREADSPLIT requires at least <fasta>:<maxErrorRate>");
		}

		// Parse optional integer parameters right-to-left, then maxErrorRate.
		// Syntax (left to right): <path>:<maxErrorRate>[:<minOverlap>[:<minFragmentLength>]]
		int   tempMinOverlap = DEFAULT_MIN_OVERLAP;
		int   tempMinFragLen = DEFAULT_MIN_FRAGMENT_LENGTH;
		float tempErrorRate;
		int   fileEndIdx;

		// Collect trailing integers (rightmost first).
		List<Integer> rightInts = new ArrayList<>();
		int scanIdx = parts.length - 1;
		while (scanIdx >= 0) {
			try {
				rightInts.add(Integer.parseInt(parts[scanIdx]));
				scanIdx--;
			} catch (NumberFormatException e) {
				break;
			}
		}

		if (scanIdx < 0) {
			throw new IllegalArgumentException(
					"LONGREADSPLIT: could not locate maxErrorRate in arguments: " + args);
		}
		try {
			tempErrorRate = Float.parseFloat(parts[scanIdx]);
		} catch (NumberFormatException e) {
			throw new IllegalArgumentException(
					"LONGREADSPLIT: invalid maxErrorRate '" + parts[scanIdx] + "'");
		}
		fileEndIdx = scanIdx;

		// rightInts (rightmost first): [minFragLen?, minOverlap?] or [minOverlap?]
		switch (rightInts.size()) {
			case 0:
				break;
			case 1:
				tempMinOverlap = rightInts.get(0);
				break;
			default:
				// rightInts.get(0) is the RIGHTMOST integer = minFragmentLength
				// rightInts.get(1) is the next     integer = minOverlap
				tempMinFragLen = rightInts.get(0);
				tempMinOverlap = rightInts.get(1);
				break;
		}

		// Reconstruct file path from parts[0..fileEndIdx-1].
		StringBuilder pathBuilder = new StringBuilder();
		for (int i = 0; i < fileEndIdx; i++) {
			if (i > 0) pathBuilder.append(':');
			pathBuilder.append(parts[i]);
		}
		if (pathBuilder.length() == 0) {
			throw new IllegalArgumentException("LONGREADSPLIT: missing FASTA file path");
		}

		maxErrorRate      = tempErrorRate;
		minOverlap        = tempMinOverlap;
		minFragmentLength = tempMinFragLen;

		adapters  = new ArrayList<>();
		kmerIndex = new HashMap<>();
		loadAdapters(new File(pathBuilder.toString()));

		if (adapters.isEmpty()) {
			throw new IllegalArgumentException(
					"LONGREADSPLIT: no adapter sequences found in " + pathBuilder);
		}
	}

	// ------------------------------------------------------------------
	// Adapter loading and k-mer indexing

	private void loadAdapters(File fastaFile) throws IOException {
		try (BufferedReader br = new BufferedReader(new FileReader(fastaFile))) {
			String line;
			StringBuilder seq = new StringBuilder();
			while ((line = br.readLine()) != null) {
				line = line.trim();
				if (line.startsWith(">")) {
					if (seq.length() > 0) {
						addAdapter(seq.toString());
						seq.setLength(0);
					}
				} else if (!line.isEmpty()) {
					seq.append(line);
				}
			}
			if (seq.length() > 0) addAdapter(seq.toString());
		}
	}

	private void addAdapter(String seq) {
		String upper = seq.toUpperCase();
		indexAdapter(upper);
		String rc = reverseComplement(upper);
		if (!rc.equals(upper)) indexAdapter(rc);
	}

	private void indexAdapter(String seq) {
		int ai = adapters.size();
		adapters.add(seq);
		for (int pos = 0; pos <= seq.length() - KMER_SIZE; pos++) {
			String kmer = seq.substring(pos, pos + KMER_SIZE);
			kmerIndex.computeIfAbsent(kmer, k -> new ArrayList<>())
			         .add(new int[]{ai, pos});
		}
	}

	// ------------------------------------------------------------------
	// Sequence utilities

	private static String reverseComplement(String seq) {
		int len = seq.length();
		StringBuilder sb = new StringBuilder(len);
		for (int i = len - 1; i >= 0; i--) sb.append(complement(seq.charAt(i)));
		return sb.toString();
	}

	private static char complement(char c) {
		return switch (c) {
			case 'A' -> 'T';
			case 'T' -> 'A';
			case 'G' -> 'C';
			case 'C' -> 'G';
			default  -> 'N';
		};
	}

	// ------------------------------------------------------------------
	// Trimmer interface

	@Override
	public FastqRecord[] processRecords(FastqRecord[] in) {
		if (in.length > 1) {
			throw new RuntimeException(
				"LONGREADSPLIT cannot be used in paired-end mode. " +
				"Long-read sequencing platforms (ONT, PacBio) produce single-end data only. " +
				"Please run Trimmomatic in SE mode.");
		}

		if (in[0] == null) return in;

		FastqRecord rec     = in[0];
		String      seq     = rec.getSequence();
		int         readLen = rec.getLength();

		// Read too short to contain any internal adapter plus two valid fragments.
		if (readLen < 2 * minOverlap + 2 * minFragmentLength) return in;

		List<int[]> hits = findInternalHits(seq, readLen);
		if (hits.isEmpty()) return in;

		// Build output fragments.
		List<FastqRecord> fragments = new ArrayList<>();
		int start   = 0;
		int total   = hits.size() + 1;
		int fragNum = 1;

		for (int[] hit : hits) {
			int hitStart = hit[0];
			int hitEnd   = hit[1];
			int fragLen  = hitStart - start;
			if (fragLen >= minFragmentLength) {
				String name = rec.getName() + "/split" + fragNum + "of" + total;
				fragments.add(new FastqRecord(rec, start, fragLen, name));
			}
			fragNum++;
			start = hitEnd;
		}

		// Last fragment.
		int lastLen = readLen - start;
		if (lastLen >= minFragmentLength) {
			String name = rec.getName() + "/split" + fragNum + "of" + total;
			fragments.add(new FastqRecord(rec, start, lastLen, name));
		}

		if (fragments.isEmpty())   return new FastqRecord[]{null};
		if (fragments.size() == 1) return new FastqRecord[]{fragments.get(0)};
		return fragments.toArray(new FastqRecord[0]);
	}

	// ------------------------------------------------------------------
	// Internal hit detection

	/**
	 * Returns a sorted, non-overlapping list of {hitStart, hitEnd} intervals
	 * for adapter matches found in the interior of {@code seq}, excluding the
	 * terminal zones at both ends.
	 */
	private List<int[]> findInternalHits(String seq, int readLen) {
		int terminalZone = minOverlap;

		// Step 1 – k-mer seeding: identify candidate adapter-start positions.
		Set<Integer> candidateStarts = new HashSet<>();
		int scanEnd = readLen - KMER_SIZE - terminalZone;
		for (int rp = terminalZone; rp <= scanEnd; rp++) {
			String kmer = seq.substring(rp, rp + KMER_SIZE);
			List<int[]> kmerHits = kmerIndex.get(kmer);
			if (kmerHits == null) continue;
			for (int[] kh : kmerHits) {
				int posInAdapter = kh[1];
				int readStart    = rp - posInAdapter;
				if (readStart >= terminalZone
						&& readStart <= readLen - terminalZone - minOverlap) {
					candidateStarts.add(readStart);
				}
			}
		}

		if (candidateStarts.isEmpty()) return Collections.emptyList();

		// Step 2 – Hamming verification at each candidate position.
		List<int[]> verified = new ArrayList<>();
		for (int readStart : candidateStarts) {
			int[] best = bestMatchAt(seq, readLen, readStart, terminalZone);
			if (best != null) verified.add(best);
		}

		if (verified.isEmpty()) return Collections.emptyList();

		// Step 3 – sort and merge overlapping intervals.
		verified.sort((a, b) -> a[0] - b[0]);
		return mergeOverlapping(verified);
	}

	/**
	 * Checks all adapter sequences at the given read-start position and returns
	 * the {hitStart, hitEnd} of the longest valid (within error budget) match,
	 * or {@code null} if no adapter matches.
	 */
	private int[] bestMatchAt(String seq, int readLen, int readStart, int terminalZone) {
		int[] best = null;
		for (int ai = 0; ai < adapters.size(); ai++) {
			String adapter    = adapters.get(ai);
			int    adapterLen = adapter.length();
			int    maxOverlap = Math.min(adapterLen, readLen - readStart - terminalZone);

			if (maxOverlap < minOverlap) continue;

			for (int overlap = maxOverlap; overlap >= minOverlap; overlap--) {
				int allowedMismatches = (int) (overlap * maxErrorRate);
				int mismatches        = 0;
				boolean exceeded      = false;

				for (int i = 0; i < overlap; i++) {
					char rc = seq.charAt(readStart + i);
					char ac = adapter.charAt(i);
					if (rc != ac && ac != 'N' && rc != 'N') {
						if (++mismatches > allowedMismatches) {
							exceeded = true;
							break;
						}
					}
				}
				if (!exceeded) {
					int hitEnd = readStart + overlap;
					if (best == null || overlap > (best[1] - best[0])) {
						best = new int[]{readStart, hitEnd};
					}
					break; // longest match per adapter found
				}
			}
		}
		return best;
	}

	/**
	 * Merge overlapping or adjacent hit intervals, preserving sorted order.
	 */
	private List<int[]> mergeOverlapping(List<int[]> sorted) {
		List<int[]> merged  = new ArrayList<>();
		int[]       current = sorted.get(0).clone();
		for (int i = 1; i < sorted.size(); i++) {
			int[] next = sorted.get(i);
			if (next[0] < current[1]) {
				current[1] = Math.max(current[1], next[1]);
			} else {
				merged.add(current);
				current = next.clone();
			}
		}
		merged.add(current);
		return merged;
	}
}
