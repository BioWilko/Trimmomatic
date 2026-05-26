package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertInstanceOf;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.mock;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.UUID;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.util.Logger;

public class TrimmerFactoryTest {

	@Test
	public void testMakeLeadingTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("LEADING:3");
		assertNotNull(trimmer);
		assertInstanceOf(LeadingTrimmer.class, trimmer);
	}

	@Test
	public void testMakeSlidingWindowTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("SLIDINGWINDOW:4:15");
		assertNotNull(trimmer);
		assertInstanceOf(SlidingWindowTrimmer.class, trimmer);
	}

	@Test
	public void testMakeTrailingTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("TRAILING:3");
		assertNotNull(trimmer);
		assertInstanceOf(TrailingTrimmer.class, trimmer);
	}

	@Test
	public void testMakeCropTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("CROP:50");
		assertNotNull(trimmer);
		assertInstanceOf(CropTrimmer.class, trimmer);
	}

	@Test
	public void testMakeHeadCropTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("HEADCROP:10");
		assertNotNull(trimmer);
		assertInstanceOf(HeadCropTrimmer.class, trimmer);
	}

	@Test
	public void testMakeTailCropTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("TAILCROP:10");
		assertNotNull(trimmer);
		assertInstanceOf(TailCropTrimmer.class, trimmer);
	}

	@Test
	public void testMakeMinLenTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("MINLEN:36");
		assertNotNull(trimmer);
		assertInstanceOf(MinLenTrimmer.class, trimmer);
	}

	@Test
	public void testMakeAvgQualTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("AVGQUAL:20");
		assertNotNull(trimmer);
		assertInstanceOf(AvgQualTrimmer.class, trimmer);
	}

	@Test
	public void testMakeMaxInfoTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("MAXINFO:50:0.5");
		assertNotNull(trimmer);
		assertInstanceOf(MaximumInformationTrimmer.class, trimmer);
	}

	@Test
	public void testMakeBaseCountTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("BASECOUNT:ACGT:5:100");
		assertNotNull(trimmer);
		assertInstanceOf(BaseCountTrimmer.class, trimmer);
	}

	@Test
	public void testMakeMaxLenTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("MAXLEN:50");
		assertNotNull(trimmer);
		assertInstanceOf(MaxLenTrimmer.class, trimmer);
	}

	@Test
	public void testMakeToPhred33Trimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("TOPHRED33");
		assertNotNull(trimmer);
		assertInstanceOf(ToPhred33Trimmer.class, trimmer);
	}

	@Test
	public void testMakeToPhred64Trimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("TOPHRED64");
		assertNotNull(trimmer);
		assertInstanceOf(ToPhred64Trimmer.class, trimmer);
	}

	@Test
	public void testMakeIlluminaClipTrimmer() throws IOException {
		File adapterFile = new File("adapters_factory_" + UUID.randomUUID() + ".fa");
		try (FileWriter writer = new FileWriter(adapterFile)) {
			writer.write(">seq\nAGATCGGAAGAGC\n");
		}

		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		String args = "ILLUMINACLIP:" + adapterFile.getName() + ":2:30:10:8:true";

		try {
			Trimmer trimmer = factory.makeTrimmer(args);
			assertNotNull(trimmer);
			assertInstanceOf(IlluminaClippingTrimmer.class, trimmer);
		} finally {
			if (!adapterFile.delete()) {
				System.gc();
				if (!adapterFile.delete())
					adapterFile.deleteOnExit();
			}
		}
	}

	@Test
	public void testMakeIlluminaClipTrimmerOptionalArgs() throws IOException {
		File adapterFile = new File("adapters_factory_opt_" + UUID.randomUUID() + ".fa");
		try (FileWriter writer = new FileWriter(adapterFile)) {
			writer.write(">seq\nAGATCGGAAGAGC\n");
		}

		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		// Args: file:seed:palindrome:simple (omitting minPrefix and keepBoth)
		String args = "ILLUMINACLIP:" + adapterFile.getName() + ":2:30:10";

		try {
			Trimmer trimmer = factory.makeTrimmer(args);
			assertNotNull(trimmer);
			assertInstanceOf(IlluminaClippingTrimmer.class, trimmer);
		} finally {
			if (!adapterFile.delete()) {
				System.gc();
				if (!adapterFile.delete())
					adapterFile.deleteOnExit();
			}
		}
	}

	@Test
	public void testMakeMaxAmbigTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("MAXAMBIG:0.2");
		assertNotNull(trimmer);
		assertInstanceOf(MaxAmbigTrimmer.class, trimmer);
	}

	@Test
	public void testMakePolyXTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("POLYX:A:10");
		assertNotNull(trimmer);
		assertInstanceOf(PolyXTrimmer.class, trimmer);
	}

	@Test
	public void testMakeLowComplexityTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("LOWCOMPLEXITY:1.0");
		assertNotNull(trimmer);
		assertInstanceOf(LowComplexityTrimmer.class, trimmer);
	}

	@Test
	public void testMakeUmiExtractTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("UMIEXTRACT:12");
		assertNotNull(trimmer);
		assertInstanceOf(UmiExtractTrimmer.class, trimmer);
	}

	@Test
	public void testUnknownTrimmer() {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		assertThrows(RuntimeException.class, () -> factory.makeTrimmer("UNKNOWN:123"));
	}

	@Test
	public void testMakeLeadingTrimmerInvalidArgs() {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		// Non-integer argument
		assertThrows(NumberFormatException.class, () -> factory.makeTrimmer("LEADING:abc"));
	}

	@Test
	public void testMakeSlidingWindowTrimmerInvalidArgs() {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		// Missing second argument
		assertThrows(ArrayIndexOutOfBoundsException.class, () -> factory.makeTrimmer("SLIDINGWINDOW:4"));
	}
}
