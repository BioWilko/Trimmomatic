package org.usadellab.trimmomatic.trim;

import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.mockito.Mockito.mock;

import org.junit.jupiter.api.Test;
import org.usadellab.trimmomatic.util.Logger;

public class TrimmerFactoryWhitespaceTest {

	@Test
	public void testLeadingWhitespace() {
		TrimmerFactory factory = new TrimmerFactory(mock(Logger.class));
		
		// " LEADING:3" should fail because " LEADING" != "LEADING"
		assertThrows(RuntimeException.class, () -> {
			factory.makeTrimmer(" LEADING:3");
		});
	}
}