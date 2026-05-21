package org.usadellab.trimmomatic.util.compression;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.EnumSet;
import java.util.Set;
import java.util.zip.GZIPOutputStream;
import java.util.zip.ZipInputStream;

import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;
import org.apache.commons.compress.compressors.bzip2.BZip2CompressorOutputStream;
import org.usadellab.trimmomatic.util.Logger;

public enum CompressionFormat {
	GZIP(".gz"), BZIP2(".bz2"), ZIP(".zip"), NONE("");

	private final String ext;
	private static Set<CompressionFormat> extSet;

	static {
		// Optimization: Use EnumSet for better performance
		extSet = EnumSet.allOf(CompressionFormat.class);
		extSet.remove(NONE);
	}

	private CompressionFormat(String ext) {
		this.ext = ext;
	}

	public static CompressionFormat forName(String name) {
		// Optimization: Avoid allocating a new String with toLowerCase()
		for (CompressionFormat cf : extSet) {
			// Case-insensitive suffix check
			if (name.regionMatches(true, name.length() - cf.ext.length(), cf.ext, 0, cf.ext.length()))
				return cf;
		}

		return NONE;
	}

	public static InputStream wrapStreamForParsing(InputStream is, String name) throws IOException {
		CompressionFormat cf = CompressionFormat.forName(name);

		return switch (cf) {
			case GZIP -> new ConcatGZIPInputStream(is);
			case BZIP2 -> new BZip2CompressorInputStream(is, true);
			case ZIP -> new ZipInputStream(is);
			case NONE -> is;
		};
	}

	public static OutputStream wrapStreamForSerializing(OutputStream os, String name, Integer compressLevel)
			throws IOException {
		CompressionFormat cf = CompressionFormat.forName(name);

		return switch (cf) {
			case GZIP -> compressLevel == null ? new GZIPOutputStream(os) : new TunableGZIPOutputStream(os, compressLevel);
			case BZIP2 -> compressLevel == null ? new BZip2CompressorOutputStream(os) : new BZip2CompressorOutputStream(os, compressLevel);
			default -> os;
		};
	}

	public static ParallelCompressor parallelCompressorForSerializing(Logger logger, String name,
			Integer compressLevel) {
		CompressionFormat cf = CompressionFormat.forName(name);

		return switch (cf) {
			case GZIP -> new GzipParallelCompressor(logger, compressLevel);
			case BZIP2 -> new Bzip2ParallelCompressor(compressLevel);
			default -> null;
		};

	}

}
