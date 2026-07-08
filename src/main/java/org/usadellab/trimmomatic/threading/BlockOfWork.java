package org.usadellab.trimmomatic.threading;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

import org.usadellab.trimmomatic.TrimStats;
import org.usadellab.trimmomatic.fastq.FastqRecord;
import org.usadellab.trimmomatic.threading.serializer.SerializedBlock;
import org.usadellab.trimmomatic.threading.serializer.Serializer;
import org.usadellab.trimmomatic.threading.trimlog.TrimLogRecord;
import org.usadellab.trimmomatic.trim.Trimmer;
import org.usadellab.trimmomatic.util.Logger;
import org.usadellab.trimmomatic.util.compression.UncompressedBlockData;

public class BlockOfWork implements Callable<BlockOfRecords> {
	private Logger logger;

	private Trimmer trimmers[];
	private BlockOfRecords bor;

	private boolean pe;
	// 0 = normal symmetric PE; 1 = R1 is the technical read; 2 = R2 is the technical read.
	private int technicalRead;
	private boolean trimLog;

	private List<SerializedBlock> blocks;
	private List<Serializer> serializers;

	private ExceptionHolder exceptionHolder;

	public BlockOfWork(Logger logger, Trimmer trimmers[], BlockOfRecords bor, boolean last, boolean pe,
			int technicalRead, boolean trimLog, List<Serializer> serializers, ExceptionHolder exceptionHolder) {
		this.logger = logger;

		this.trimmers = trimmers;
		this.bor = bor;

		this.pe = pe;
		this.technicalRead = technicalRead;
		this.trimLog = trimLog;

		blocks = new ArrayList<SerializedBlock>();

		for (@SuppressWarnings("unused")
		Serializer serializer : serializers) {
			blocks.add(new SerializedBlock(last));
		}

		this.serializers = serializers;

		this.exceptionHolder = exceptionHolder;
	}

	public List<SerializedBlock> getBlocks() {
		return blocks;
	}

	private TrimLogRecord makeTrimLogRec(FastqRecord rec, FastqRecord originalRec) {
		int length = 0;
		int startPos = 0;
		int endPos = 0;
		int trimTail = 0;

		if (rec != null) {
			length = rec.getLength();
			startPos = rec.getHeadPos();
			endPos = length + startPos;
			trimTail = originalRec.getLength() - endPos;
		}

		return new TrimLogRecord(originalRec.getName(), length, startPos, endPos, trimTail);
	}

	public void processPE() throws Exception {
		List<FastqRecord> originalRecs1 = bor.takeOriginalRecs1();
		List<FastqRecord> originalRecs2 = bor.takeOriginalRecs2();

		int len1 = originalRecs1.size();
		int len2 = originalRecs2.size();

		if (len1 == 0 && len2 == 0) {
			for (SerializedBlock block : blocks)
				block.setUncompressedData(null);

			for (Serializer ser : serializers)
				ser.pollCompressible();

			bor.setTrimLogRecs(null);
			return;
		}

		int len = len1 < len2 ? len1 : len2;

		FastqRecord originalRecs[] = new FastqRecord[2];

		List<FastqRecord> trimmedRecs1P = new ArrayList<FastqRecord>();
		List<FastqRecord> trimmedRecs1U = new ArrayList<FastqRecord>();
		List<FastqRecord> trimmedRecs2P = new ArrayList<FastqRecord>();
		List<FastqRecord> trimmedRecs2U = new ArrayList<FastqRecord>();

		List<TrimLogRecord> trimLogList = null;
		if (trimLog)
			trimLogList = new ArrayList<TrimLogRecord>();

		TrimStats stats = new TrimStats();

		for (int i = 0; i < len; i++) {
			originalRecs[0] = originalRecs1.get(i);
			originalRecs[1] = originalRecs2.get(i);

			if (technicalRead != 0) {
				// Asymmetric PE: one read is a technical read (barcode/UMI) that passes
				// through completely untouched; all trimmers run on the biological read only.
				// If the biological read is dropped, both reads are discarded — the technical
				// read NEVER goes to the unpaired output.
				int bioIdx  = (technicalRead == 1) ? 1 : 0;
				int techIdx = (technicalRead == 1) ? 0 : 1;

				// Run all trimmers on the biological read in SE mode (single-element array).
				FastqRecord[] bioArray = { originalRecs[bioIdx] };
				for (int j = 0; j < trimmers.length; j++) {
					try {
						bioArray = trimmers[j].processRecords(bioArray);
					} catch (RuntimeException e) {
						logger.errorln("Exception processing bio read: " + originalRecs[bioIdx].getName());
						throw e;
					}
				}

				FastqRecord bioResult  = bioArray[0];
				FastqRecord techRecord = originalRecs[techIdx];

				if (bioResult != null) {
					// Bio read survived → emit pair; tech read is passed through verbatim.
					if (technicalRead == 1) {
						trimmedRecs1P.add(techRecord); // R1 = technical
						trimmedRecs2P.add(bioResult);  // R2 = biological
					} else {
						trimmedRecs1P.add(bioResult);  // R1 = biological
						trimmedRecs2P.add(techRecord); // R2 = technical
					}
				}
				// Bio read dropped → both discarded. Technical read NEVER goes to unpaired.

				// Stats: pair survives iff bio survives; they always move together.
				FastqRecord[] recsForStats = new FastqRecord[2];
				recsForStats[bioIdx]  = bioResult;
				recsForStats[techIdx] = (bioResult != null) ? techRecord : null;
				stats.logPair(originalRecs, recsForStats);

				if (trimLog) {
					trimLogList.add(makeTrimLogRec(recsForStats[0], originalRecs[0]));
					trimLogList.add(makeTrimLogRec(recsForStats[1], originalRecs[1]));
				}
			} else {
				// Normal symmetric PE processing.
				FastqRecord recs[] = originalRecs;

				for (int j = 0; j < trimmers.length; j++) {
					try {
						recs = trimmers[j].processRecords(recs);
					} catch (RuntimeException e) {
						logger.errorln("Exception processing reads: " + originalRecs[0].getName() + " and "
								+ originalRecs[1].getName());
						throw e;
					}
				}

				if (recs[0] != null && recs[1] != null) {
					trimmedRecs1P.add(recs[0]);
					trimmedRecs2P.add(recs[1]);
				} else if (recs[0] != null)
					trimmedRecs1U.add(recs[0]);
				else if (recs[1] != null)
					trimmedRecs2U.add(recs[1]);

				stats.logPair(originalRecs, recs);

				if (trimLog) {
					if (originalRecs[0] != null)
						trimLogList.add(makeTrimLogRec(recs[0], originalRecs[0]));
					if (originalRecs[1] != null)
						trimLogList.add(makeTrimLogRec(recs[1], originalRecs[1]));
				}
			}
		}

		blocks.get(0).setUncompressedData(new UncompressedBlockData(trimmedRecs1P));
		blocks.get(1).setUncompressedData(new UncompressedBlockData(trimmedRecs1U));
		blocks.get(2).setUncompressedData(new UncompressedBlockData(trimmedRecs2P));
		blocks.get(3).setUncompressedData(new UncompressedBlockData(trimmedRecs2U));

		for (Serializer ser : serializers)
			ser.pollCompressible();

		bor.setTrimLogRecs(trimLogList);
		bor.setStats(stats);

	}

	public void processSE() throws Exception {
		List<FastqRecord> originalRecsL = bor.takeOriginalRecs1();

		int len = originalRecsL.size();

		if (len == 0) {
			blocks.get(0).setUncompressedData(null);
			serializers.get(0).pollCompressible();

			bor.setTrimLogRecs(null);
			return;
		}

		FastqRecord originalRecs[] = new FastqRecord[1];

		List<FastqRecord> trimmedRecs = new ArrayList<FastqRecord>();

		List<TrimLogRecord> trimLogList = null;
		if (trimLog)
			trimLogList = new ArrayList<TrimLogRecord>();

		TrimStats stats = new TrimStats();

		for (int i = 0; i < len; i++) {
			originalRecs[0] = originalRecsL.get(i);
			FastqRecord recs[] = originalRecs;

			for (int j = 0; j < trimmers.length; j++) {
				try {
					recs = trimmers[j].processRecords(recs);
				} catch (RuntimeException e) {
					logger.errorln("Exception processing read: " + originalRecs[0].getName());
					e.printStackTrace();
					throw e;
				}
			}

			for (FastqRecord rec : recs) {
				if (rec != null)
					trimmedRecs.add(rec);
			}

			stats.logPair(originalRecs, recs);

			if (trimLog && originalRecs[0] != null) {
				for (FastqRecord rec : recs) {
					trimLogList.add(makeTrimLogRec(rec, originalRecs[0]));
				}
			}
		}

		blocks.get(0).setUncompressedData(new UncompressedBlockData(trimmedRecs));
		serializers.get(0).pollCompressible();

		bor.setTrimLogRecs(trimLogList);
		bor.setStats(stats);
	}

	public BlockOfRecords process() throws Exception {
		if (pe)
			processPE();
		else
			processSE();

		return bor;
	}

	@Override
	public BlockOfRecords call() {
		try {
			process();
			// throw new Exception("FakeBreak");
		} catch (Exception e) {
			Exception pe = new Exception("Pipeline Exception", e);
			exceptionHolder.setException(pe);
		}

		return bor;
	}

}
