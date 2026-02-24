package org.usadellab.trimmomatic.fastq;

import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.mockito.ArgumentMatchers.contains;
import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.verify;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.util.Logger;

public class PairingValidatorTest {

	private PairingValidator validator;
	private Logger logger;

	@BeforeEach
	public void setUp() {
		logger = mock(Logger.class);
		validator = new PairingValidator(logger);
	}

	private FastqRecord makeRecord(String name) {
		// 150bp dummy record
		return new FastqRecord(name, "ACGT", "", "IIII", 33);
	}

	@Test
	public void testValidateIdenticalNames() {
		FastqRecord r1 = makeRecord("seq1");
		FastqRecord r2 = makeRecord("seq1");
		assertTrue(validator.validatePair(r1, r2));
	}

	@Test
	public void testValidateCasavaSuffix() {
		// Old Casava /1 and /2 style
		FastqRecord r1 = makeRecord("seq1/1");
		FastqRecord r2 = makeRecord("seq1/2");
		assertTrue(validator.validatePair(r1, r2));
	}

	@Test
	public void testValidateCasavaSpace() {
		// Casava 1.8+ space style
		FastqRecord r1 = makeRecord("seq1 1:Y:0:ATCACG");
		FastqRecord r2 = makeRecord("seq1 2:Y:0:ATCACG");
		assertTrue(validator.validatePair(r1, r2));
	}

	@Test
	public void testValidateMismatchNames() {
		FastqRecord r1 = makeRecord("seqA");
		FastqRecord r2 = makeRecord("seqB");
		assertFalse(validator.validatePair(r1, r2));
		verify(logger).warnln(contains("Pair validation failed"));
	}

	@Test
	public void testValidateMismatchSuffix() {
		FastqRecord r1 = makeRecord("seq1/1");
		FastqRecord r2 = makeRecord("seq1/3"); // Invalid suffix
		assertFalse(validator.validatePair(r1, r2));
	}

	@Test
	public void testValidateNullRecord() {
		FastqRecord r1 = makeRecord("seq1");
		// Missing reverse read
		assertFalse(validator.validatePair(r1, null));
		verify(logger).warnln(contains("No more reverse reads"));
		
		// Missing forward read
		assertFalse(validator.validatePair(null, r1));
	}
}