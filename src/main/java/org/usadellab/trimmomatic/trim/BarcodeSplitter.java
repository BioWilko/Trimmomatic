/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.usadellab.trimmomatic.trim;

import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Map;

import org.usadellab.trimmomatic.fastq.FastqRecord;

/**
 *
 * @author marc
 */
public class BarcodeSplitter extends AbstractSingleRecordTrimmer {

	private HashMap<String, String> barcodes;
	private Barcode[] barcodeArray;
	private int maxMisMatch = 0;
	private boolean clipOffBarcodes = true;

	public BarcodeSplitter(HashMap<String, String> barcodes, int mism, boolean clip) {
		this.barcodes = barcodes;
		this.maxMisMatch = mism;
		this.clipOffBarcodes = clip;

		// Optimization: Flatten map to array to avoid Iterator allocation per record
		this.barcodeArray = new Barcode[barcodes.size()];
		int i = 0;
		for (Map.Entry<String, String> entry : barcodes.entrySet()) {
			this.barcodeArray[i++] = new Barcode(entry.getKey(), entry.getValue());
		}

		// Optimization: Sort barcodes by length (descending) to ensure longest match is found first,
		// and by label for determinism.
		Arrays.sort(this.barcodeArray, Comparator.comparingInt((Barcode b) -> b.length).reversed().thenComparing(b -> b.label));
	}

	public HashMap<String, String> getBarcodeMap() {
		return barcodes;
	}

	@Override
	public FastqRecord processRecord(FastqRecord entry) {

		String seq = entry.getSequence();
		int seqLen = seq.length();

		// Optimization: Unswitch loop to avoid checking maxMisMatch inside the loop
		if (maxMisMatch == 0) {
			for (Barcode b : barcodeArray) {
				if (seqLen < b.length) continue;

				if (seq.startsWith(b.sequence)) {
					return createOutputRecord(entry, b.label, b.length, seqLen);
				}
			}
		} else {
			for (Barcode b : barcodeArray) {
				if (seqLen < b.length) continue;

				if (getMisMatches(b.sequenceChars, seq, maxMisMatch) <= maxMisMatch) {
					return createOutputRecord(entry, b.label, b.length, seqLen);
				}
			}
		}
		// entry did not match any barcodes - write to UNKNOWN
		entry.setBarcodeLabel("UNKNOWN");
		return entry;
	}

	private FastqRecord createOutputRecord(FastqRecord entry, String label, int codeLength, int seqLength) {
		if (clipOffBarcodes) {
			FastqRecord outentry = new FastqRecord(entry, codeLength, seqLength - codeLength);
			outentry.setBarcodeLabel(label);
			return outentry;
		} else {
			entry.setBarcodeLabel(label);
			return entry;
		}
	}

	private static int getMisMatches(char[] barcodeChars, String b, int maxAllowed) {
		int mismatch = 0;
		for (int i = 0; i < barcodeChars.length; i++) {
			if (barcodeChars[i] != b.charAt(i)) {
				mismatch++;
				if (mismatch > maxAllowed) return mismatch;
			}
		}
		return mismatch;
	}

	private static class Barcode {
		final String label;
		final String sequence;
		final char[] sequenceChars;
		final int length;

		Barcode(String label, String sequence) {
			this.label = label;
			this.sequence = sequence;
			this.sequenceChars = sequence.toCharArray();
			this.length = sequence.length();
		}
	}
}
