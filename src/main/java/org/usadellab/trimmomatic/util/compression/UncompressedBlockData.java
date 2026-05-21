package org.usadellab.trimmomatic.util.compression;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.List;

import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.fastq.FastqSerializer;

public final class UncompressedBlockData implements BlockData {
	private byte[] data;
	private int dataLength;

	// Exposes the internal buffer so we can avoid the extra copy that
	// ByteArrayOutputStream.toByteArray() would otherwise make.
	private static final class ExposedByteArrayOutputStream extends ByteArrayOutputStream {
		ExposedByteArrayOutputStream(int capacity) {
			super(capacity);
		}

		byte[] getBuffer() {
			return buf;
		}

		int getCount() {
			return count;
		}
	}

	public UncompressedBlockData(List<FastqRecord> recs) {
		int estimated = 0;
		for (FastqRecord rec : recs)
			estimated += rec.getRecordLength();

		ExposedByteArrayOutputStream stream = new ExposedByteArrayOutputStream(estimated);
		FastqSerializer ser = new FastqSerializer();

		try {
			ser.open(stream);
			for (FastqRecord rec : recs)
				ser.writeRecord(rec);
			ser.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		// If the buffer was exactly filled (common when estimate is accurate) reuse it
		// directly; otherwise trim to the actual written length.
		int count = stream.getCount();
		byte[] buf = stream.getBuffer();
		if (count == buf.length) {
			data = buf;
		} else {
			data = new byte[count];
			System.arraycopy(buf, 0, data, 0, count);
		}
		dataLength = count;
	}

	public byte[] data() {
		return data;
	}
}
