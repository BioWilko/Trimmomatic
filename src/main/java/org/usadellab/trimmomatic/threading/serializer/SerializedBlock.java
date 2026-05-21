package org.usadellab.trimmomatic.threading.serializer;

import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;

import org.usadellab.trimmomatic.util.compression.BlockData;
import org.usadellab.trimmomatic.util.compression.UncompressedBlockData;

public class SerializedBlock {
    boolean last;

    private final CompletableFuture<UncompressedBlockData> ucFuture = new CompletableFuture<>();
    private final CompletableFuture<BlockData> dataFuture = new CompletableFuture<>();

    public SerializedBlock(boolean last) {
        this.last = last;
    }

    public boolean isLast() {
        return last;
    }

    public void setUncompressedData(UncompressedBlockData ucData) {
        ucFuture.complete(ucData);
    }

    public UncompressedBlockData getUncompressedData() {
        return ucFuture.getNow(null);
    }

    public boolean isCompressible() {
        return ucFuture.isDone();
    }

    public void setData(BlockData data) {
        dataFuture.complete(data);
    }

    public boolean isDone() {
        return dataFuture.isDone();
    }

    public BlockData get() throws InterruptedException, ExecutionException {
        return dataFuture.get();
    }

}
