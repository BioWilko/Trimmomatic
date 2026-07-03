package org.usadellab.trimmomatic.util.compression;

public record GzipBlockData(byte[] data) implements BlockData {
	public GzipBlockData {
		if (data == null)
			throw new NullPointerException("GzipBlockData: cannot be null");
	}
}
