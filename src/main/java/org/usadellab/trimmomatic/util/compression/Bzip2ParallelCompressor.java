package org.usadellab.trimmomatic.util.compression;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.Collections;

import org.apache.commons.compress.compressors.bzip2.BZip2CompressorOutputStream;

public class Bzip2ParallelCompressor implements ParallelCompressor {
	private static final int DEFAULT_COMPRESSION_LEVEL = 9;

	private int compressLevel;

	public Bzip2ParallelCompressor(Integer compressLevel) {
		this.compressLevel = compressLevel == null ? DEFAULT_COMPRESSION_LEVEL : compressLevel;
	}

	@Override
	public BlockOutputStream wrapAndWriteHeader(OutputStream stream) throws IOException {
		return new BlockOutputStream(stream);
	}

	@Override
	public Bzip2BlockData compress(UncompressedBlockData previous, UncompressedBlockData current) throws Exception {
		if (current == null)
			return null;

		byte[] rawData = current.data();

		if (rawData.length == 0)
			return new Bzip2BlockData(new byte[0], 0, Collections.emptyList());

		try {
			// Optimization: Pre-allocate buffer to avoid resizing churn.
			ByteArrayOutputStream baos = new ByteArrayOutputStream(rawData.length / 2);
			try (BZip2CompressorOutputStream bzos = new BZip2CompressorOutputStream(baos, compressLevel)) {
				bzos.write(rawData);
			}

			byte[] out = baos.toByteArray();

			return new Bzip2BlockData(out, out.length * 8L, Collections.emptyList());
		} catch (IOException e) {
			throw new Exception("Bzip2 compression exception", e);
		}
	}
}
