package org.usadellab.trimmomatic.util.compression;

import java.util.List;

public record Bzip2BlockData(byte[] data, long bitCount, List<Integer> blockCRCs) implements BlockData {
	public Bzip2BlockData {
		if (data == null)
			throw new NullPointerException("Bzip2BlockData: cannot be null");
	}
}
