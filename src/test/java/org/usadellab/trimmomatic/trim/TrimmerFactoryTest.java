package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;
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
		assertTrue(trimmer instanceof LeadingTrimmer);
	}

	@Test
	public void testMakeSlidingWindowTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("SLIDINGWINDOW:4:15");
		assertNotNull(trimmer);
		assertTrue(trimmer instanceof SlidingWindowTrimmer);
	}
	
	@Test
	public void testMakeTrailingTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("TRAILING:3");
		assertNotNull(trimmer);
		assertTrue(trimmer instanceof TrailingTrimmer);
	}

	@Test
	public void testMakeCropTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("CROP:50");
		assertNotNull(trimmer);
		assertTrue(trimmer instanceof CropTrimmer);
	}

	@Test
	public void testMakeHeadCropTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("HEADCROP:10");
		assertNotNull(trimmer);
		assertTrue(trimmer instanceof HeadCropTrimmer);
	}

	@Test
	public void testMakeTailCropTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("TAILCROP:10");
		assertNotNull(trimmer);
		assertTrue(trimmer instanceof TailCropTrimmer);
	}

	@Test
	public void testMakeMinLenTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("MINLEN:36");
		assertNotNull(trimmer);
		assertTrue(trimmer instanceof MinLenTrimmer);
	}

	@Test
	public void testMakeAvgQualTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("AVGQUAL:20");
		assertNotNull(trimmer);
		assertTrue(trimmer instanceof AvgQualTrimmer);
	}

	@Test
	public void testMakeMaxInfoTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("MAXINFO:50:0.5");
		assertNotNull(trimmer);
		assertTrue(trimmer instanceof MaximumInformationTrimmer);
	}

	@Test
	public void testMakeBaseCountTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("BASECOUNT:ACGT:5:100");
		assertNotNull(trimmer);
		assertTrue(trimmer instanceof BaseCountTrimmer);
	}

	@Test
	public void testMakeMaxLenTrimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("MAXLEN:50");
		assertNotNull(trimmer);
		assertTrue(trimmer instanceof MaxLenTrimmer);
	}

	@Test
	public void testMakeToPhred33Trimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("TOPHRED33");
		assertNotNull(trimmer);
		assertTrue(trimmer instanceof ToPhred33Trimmer);
	}

	@Test
	public void testMakeToPhred64Trimmer() throws IOException {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		Trimmer trimmer = factory.makeTrimmer("TOPHRED64");
		assertNotNull(trimmer);
		assertTrue(trimmer instanceof ToPhred64Trimmer);
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
			assertTrue(trimmer instanceof IlluminaClippingTrimmer);
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
			assertTrue(trimmer instanceof IlluminaClippingTrimmer);
		} finally {
			if (!adapterFile.delete()) {
				System.gc();
				if (!adapterFile.delete())
					adapterFile.deleteOnExit();
			}
		}
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