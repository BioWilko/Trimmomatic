package org.usadellab.trimmomatic.threading.pipeline;

import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.ThreadFactory;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import org.usadellab.trimmomatic.threading.BlockOfRecords;
import org.usadellab.trimmomatic.threading.BlockOfWork;
import org.usadellab.trimmomatic.threading.ExceptionHolder;

public class ThreadedPipeline extends Pipeline {
	private ExceptionHolder exceptionHolder;

	ArrayBlockingQueue<Runnable> taskQueue;
	ThreadPoolExecutor taskExec;

	public ThreadedPipeline(int threads, ExceptionHolder exceptionHolder) {
		this.exceptionHolder = exceptionHolder;

		taskQueue = new ArrayBlockingQueue<Runnable>(threads * 8);

		ThreadFactory tf = new ThreadFactory() {
			public Thread newThread(Runnable r) {
				Thread t = Executors.defaultThreadFactory().newThread(r);
				t.setDaemon(true);
				return t;
			}
		};

		taskExec = new ThreadPoolExecutor(threads, threads, 0, TimeUnit.SECONDS, taskQueue, tf);

		// Block the caller instead of throwing RejectedExecutionException when the
		// queue is full. This replaces the old sleep(100) busy-wait loop.
		taskExec.setRejectedExecutionHandler((r, executor) -> {
			try {
				executor.getQueue().put(r);
			} catch (InterruptedException e) {
				Thread.currentThread().interrupt();
				throw new RuntimeException(e);
			}
		});
	}

	public Future<BlockOfRecords> submit(BlockOfWork work) throws Exception {
		exceptionHolder.rethrow();
		return taskExec.submit(work);
	}

	public void close() throws InterruptedException {
		taskExec.shutdown();
		taskExec.awaitTermination(1, TimeUnit.HOURS);
	}

}
