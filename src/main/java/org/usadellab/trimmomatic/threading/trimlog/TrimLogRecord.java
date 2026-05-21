package org.usadellab.trimmomatic.threading.trimlog;

public record TrimLogRecord(String readName, int length, int startPos, int endPos, int trimTail) {}
