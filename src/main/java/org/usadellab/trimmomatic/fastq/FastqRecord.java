package org.usadellab.trimmomatic.fastq;

public class FastqRecord {
	private String name;
	private String sequence;
	private String comment;
	private String quality;
	private String barcodeLabel = null;

	private int phredOffset;
	private int headPos;

	public FastqRecord(String name, String sequence, String comment, String quality, int phredOffset) {
		this.name = name;
		this.sequence = sequence;
		this.comment = comment;
		this.quality = quality;

		this.phredOffset = phredOffset;
		headPos = 0;

		if (sequence.length() != quality.length())
			throw new RuntimeException(
					"Sequence and quality length don't match: '" + sequence + "' vs '" + quality + "'");
	}

	public FastqRecord(FastqRecord base, int headPos, int length) {
		if (headPos < 0)
			throw new RuntimeException("Attempting invalid trim on " + base.name + " with length "
					+ base.sequence.length() + ": Wanted " + headPos + " to " + (headPos + length));

		int availableLength = base.getSequence().length();
		if (headPos + length > availableLength)
			length = availableLength - headPos;

		this.sequence = base.sequence.substring(headPos, headPos + length);
		this.quality = base.quality.substring(headPos, headPos + length);
		this.name = base.name;
		this.comment = base.comment;
		this.phredOffset = base.phredOffset;
		this.headPos = base.headPos + headPos;

		this.barcodeLabel = base.barcodeLabel;
	}

	public FastqRecord(FastqRecord base, String sequence, String quality, int phredOffset) {
		this.sequence = sequence;
		this.quality = quality;
		this.name = base.name;
		this.comment = base.comment;
		this.headPos = base.headPos;
		this.phredOffset = phredOffset;

		this.barcodeLabel = base.barcodeLabel;
	}

	public static FastqRecord make(String name, String sequence, int quality) {
		StringBuilder qualityBuilder = new StringBuilder();
		for (int i = 0; i < sequence.length(); i++)
			qualityBuilder.append((char) (33 + quality));

		return new FastqRecord(name, sequence, name, qualityBuilder.toString(), 33);
	}

	public static FastqRecord make(String name, String sequence) {
		return make(name, sequence, 40);
	}

	public String getName() {
		return name;
	}

	public String getSequence() {
		return sequence;
	}

	public String getBarcodeLabel() {
		return barcodeLabel;
	}

	public void setBarcodeLabel(String barcodeLabel) {
		this.barcodeLabel = barcodeLabel;
	}

	public String getComment() {
		return comment;
	}

	public String getQuality() {
		return quality;
	}

	public int getPhredOffset() {
		return phredOffset;
	}

	void setPhredOffset(int phredOffset) {
		this.phredOffset = phredOffset;
	}

	public int getHeadPos() {
		return headPos;
	}

	public int[] getQualityAsInteger(boolean zeroNs) {
		int arr[] = new int[quality.length()];

		for (int i = 0; i < quality.length(); i++) {
			if (zeroNs && sequence.charAt(i) == 'N')
				arr[i] = 0;
			else
				arr[i] = quality.charAt(i) - phredOffset;
		}

		return arr;
	}

	private static int RECORD_ADDED_LENGTH = 10; // 4 newlines (2 characters each), @ and +

	public int getRecordLength() {
		return this.name.length() + this.sequence.length() + this.comment.length() + this.quality.length()
				+ RECORD_ADDED_LENGTH;
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((comment == null) ? 0 : comment.hashCode());
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		result = prime * result + phredOffset;
		result = prime * result + ((quality == null) ? 0 : quality.hashCode());
		result = prime * result + ((sequence == null) ? 0 : sequence.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		FastqRecord other = (FastqRecord) obj;
		if (comment == null) {
			if (other.comment != null)
				return false;
		} else if (!comment.equals(other.comment))
			return false;
		if (name == null) {
			if (other.name != null)
				return false;
		} else if (!name.equals(other.name))
			return false;
		if (phredOffset != other.phredOffset)
			return false;
		if (quality == null) {
			if (other.quality != null)
				return false;
		} else if (!quality.equals(other.quality))
			return false;
		if (sequence == null) {
			if (other.sequence != null)
				return false;
		} else if (!sequence.equals(other.sequence))
			return false;
		return true;
	}

}
