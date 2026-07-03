package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.util.Logger;

public class TrimmerFactoryEmptyTest {

	@Test
	public void testEmptyString() {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		// Empty string -> trimmerName is empty -> Unknown trimmer
		assertThrows(RuntimeException.class, () -> {
			factory.makeTrimmer("");
		});
	}
}