package org.usadellab.trimmomatic.util.compression;

public sealed interface BlockData permits GzipBlockData, Bzip2BlockData, UncompressedBlockData {
	byte[] data();
}
