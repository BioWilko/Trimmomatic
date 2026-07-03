package org.usadellab.trimmomatic.util;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

public class LoggerTest {

	private final ByteArrayOutputStream outContent = new ByteArrayOutputStream();
	private final ByteArrayOutputStream errContent = new ByteArrayOutputStream();
	private final PrintStream originalOut = System.out;
	private final PrintStream originalErr = System.err;

	@BeforeEach
	public void setUpStreams() {
		System.setOut(new PrintStream(outContent));
		System.setErr(new PrintStream(errContent));
	}

	@AfterEach
	public void restoreStreams() {
		System.setOut(originalOut);
		System.setErr(originalErr);
	}

	@Test
	public void testLogOutput() {
		Logger logger = new Logger(true, true, true); // Verbose
		logger.infoln("Info message");
		logger.warnln("Warning message");
		logger.errorln("Error message");

		// Check if output contains messages
		String errOutput = errContent.toString();
		assertTrue(errOutput.contains("Info message"), "Info message should be present in verbose mode");
		assertTrue(errOutput.contains("Warning message"), "Warning message should be present");
		assertTrue(errOutput.contains("Error message"), "Error message should be present");
	}
}