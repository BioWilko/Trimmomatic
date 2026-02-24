package org.usadellab.trimmomatic.util;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.ByteArrayInputStream;
import java.io.IOException;

import org.junit.jupiter.api.Test;

public class PositionTrackingInputStreamTest {

	@Test
	public void testProgress() throws IOException {
		byte[] data = new byte[100];
		PositionTrackingInputStream stream = new PositionTrackingInputStream(new ByteArrayInputStream(data), 100);
		
		assertEquals(0, stream.getProgressPercentage());
		
		stream.read(); // Read 1 byte
		assertEquals(1, stream.getProgressPercentage());
		
		byte[] buf = new byte[49];
		stream.read(buf); // Read 49 bytes -> Total 50
		assertEquals(50, stream.getProgressPercentage());
		
		stream.read(new byte[50]); // Read 50 bytes -> Total 100
		assertEquals(100, stream.getProgressPercentage());
		
		stream.close();
	}

	@Test
	public void testSkip() throws IOException {
		byte[] data = new byte[100];
		PositionTrackingInputStream stream = new PositionTrackingInputStream(new ByteArrayInputStream(data), 100);
		
		assertEquals(0, stream.getProgressPercentage());
		
		stream.skip(50);
		assertEquals(50, stream.getProgressPercentage());
		
		stream.skip(50);
		assertEquals(100, stream.getProgressPercentage());
		
		stream.close();
	}
	
	@Test
	public void testReadRange() throws IOException {
		byte[] data = new byte[100];
		PositionTrackingInputStream stream = new PositionTrackingInputStream(new ByteArrayInputStream(data), 100);
		
		byte[] buf = new byte[100];
		stream.read(buf, 0, 25);
		assertEquals(25, stream.getProgressPercentage());
		
		stream.close();
	}
}