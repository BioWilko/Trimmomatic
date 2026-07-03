package org.usadellab.trimmomatic.util.compression;

import java.io.IOException;
import java.io.OutputStream;

public interface ParallelCompressor {
	BlockOutputStream wrapAndWriteHeader(OutputStream stream) throws IOException;

	default void updateChecksumPreCompression(UncompressedBlockData current) {}

	BlockData compress(UncompressedBlockData previous, UncompressedBlockData current) throws Exception;

	default void updateChecksumPostCompression(BlockData current) {}

	default void writeTrailer(BlockOutputStream stream) throws IOException {}
}
