//******************************************************************************
//
// File:    WorkerLongStrideForLoop.java
// Package: edu.rit.pj
// Unit:    Class edu.rit.pj.WorkerLongStrideForLoop
//
// This Java source file is copyright (C) 2010 by Alan Kaminsky. All rights
// reserved. For further information, contact the author, Alan Kaminsky, at
// ark@cs.rit.edu.
//
// This Java source file is part of the Parallel Java Library ("PJ"). PJ is free
// software; you can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// PJ is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// Linking this library statically or dynamically with other modules is making a
// combined work based on this library. Thus, the terms and conditions of the GNU
// General Public License cover the whole combination.
//
// As a special exception, the copyright holders of this library give you
// permission to link this library with independent modules to produce an
// executable, regardless of the license terms of these independent modules, and
// to copy and distribute the resulting executable under terms of your choice,
// provided that you also meet, for each linked independent module, the terms
// and conditions of the license of that module. An independent module is a module
// which is not derived from or based on this library. If you modify this library,
// you may extend this exception to your version of the library, but you are not
// obligated to do so. If you do not wish to do so, delete this exception
// statement from your version.
//
// A copy of the GNU General Public License is provided in the file gpl.txt. You
// may also obtain a copy of the GNU General Public License on the World Wide
// Web at http://www.gnu.org/licenses/gpl.html.
//
//******************************************************************************
package edu.rit.pj;

import java.io.IOException;

import edu.rit.mp.ObjectBuf;
import edu.rit.mp.buf.ObjectItemBuf;
import edu.rit.util.LongRange;
import edu.rit.util.Range;

/**
 * Class WorkerLongStrideForLoop is the abstract base class for one variation of
 * a worker for loop that is executed inside a {@linkplain WorkerRegion}. The
 * loop index data type is <code>long</code>. The loop stride is explicitly
 * specified.
 * <P>
 * To execute a worker for loop, create a {@linkplain WorkerRegion} object;
 * create an instance of a concrete subclass of class WorkerLongStrideForLoop;
 * and pass this instance to the worker region's <code>execute()</code> method.
 * Either every worker team thread must call the worker region's
 * <code>execute()</code> method with identical arguments, or every thread must not
 * call the <code>execute()</code> method. You can do all this using an anonymous
 * inner class; for example:
 * <PRE>
 *     new WorkerRegion()
 *         {
 *         . . .
 *         public void run()
 *             {
 *             . . .
 *             execute (0L, 98L, 2L, new WorkerLongStrideForLoop()
 *                 {
 *                 // Thread local variable declarations
 *                 . . .
 *                 public void start()
 *                     {
 *                     // Per-thread pre-loop initialization code
 *                     . . .
 *                     }
 *                 public void run (long first, long last, long stride)
 *                     {
 *                     // Loop code
 *                     . . .
 *                     }
 *                 public void finish()
 *                     {
 *                     // Per-thread post-loop finalization code
 *                     . . .
 *                     }
 *                 });
 *             }
 *         . . .
 *         }
 * </PRE>
 * <P>
 * In each process of a cluster parallel program, the worker team has one or
 * more worker threads. Every worker thread in every process has a unique worker
 * tag, going from tag 0 for the first worker thread in the first process to tag
 * <I>K</I>&minus;1 for the last worker thread in the last process, where
 * <I>K</I> is the total number of worker threads in all the processes. In
 * addition, in one process there is a master thread. The worker and master
 * threads all call the worker region's <code>execute()</code> method to execute the
 * worker for loop. However, the worker and master threads differ in their
 * actions.
 * <P>
 * The master thread does the following. The master obtains the worker for
 * loop's schedule as returned by the <code>schedule()</code> method. The range of
 * loop indexes is divided into "chunks" and the chunks are apportioned among
 * the workers in accordance with the schedule. The master repeatedly sends
 * "tasks" to the workers and receives "responses" from the workers. To send a
 * task to a particular worker, the master (1) sends a message containing the
 * chunk index range to the worker's process; and (2) calls the worker for
 * loop's <code>sendTaskInput()</code> method. This method's default implementation
 * does nothing, but it can be overridden to send additional task input data to
 * the worker. To receive a response from a particular worker, the master (1)
 * receives a message containing the chunk index range from the worker's
 * process; and (2) calls the worker for loop's <code>receiveTaskOutput()</code>
 * method. This method's default implementation does nothing, but it can be
 * overridden to receive additional task output data from the worker. Once all
 * tasks have been sent to the workers and all responses have been received from
 * the workers, the master returns from the worker region's <code>execute()</code>
 * method.
 * <P>
 * Each worker thread does the following. The worker calls the worker for loop's
 * <code>start()</code> method once before beginning any loop iterations. The worker
 * repeatedly receives tasks from the master and sends responses to the master.
 * To receive a task from the master, the worker (1) receives a message
 * containing the chunk index range from the master's process; and (2) calls the
 * worker for loop's <code>receiveTaskInput()</code> method. This method's default
 * implementation does nothing, but it can be overridden to receive additional
 * task input data from the master. The worker now calls the worker for loop's
 * <code>run()</code> method, passing in the chunk index range lower and upper
 * bounds. When the <code>run()</code> method returns, the worker sends the response
 * to the master. To send the response, the worker (1) sends a message
 * containing the chunk index range to the master's process; and (2) calls the
 * worker for loop's <code>sendTaskOutput()</code> method. This method's default
 * implementation does nothing, but it can be overridden to send additional task
 * output data to the master. Once all tasks have been received from the master
 * and all responses have been sent to the master, the worker calls the worker
 * for loop's <code>finish()</code> method. (Unlike a {@linkplain ParallelTeam}'s
 * threads, the workers do <I>not</I> synchronize with each other at a barrier
 * at this point.) The worker then returns from the worker region's
 * <code>execute()</code> method.
 * <P>
 * If the worker for loop has a fixed schedule (in which there is exactly one
 * chunk with a predetermined index range for each worker), then the messages
 * containing the chunk index range are omitted, and each worker gets its chunk
 * index range directly from the fixed schedule. However, the task input data
 * (if any) and task output data (if any) are still sent and received.
 * <P>
 * Each message described above is sent with a message tag equal to
 * <I>W</I>+<I>T</I>, where <I>W</I> is the worker index and <I>T</I> is the
 * "tag offset." The tag offset is <code>Integer.MIN_VALUE</code> by default, but
 * this can be changed by overriding the <code>tagOffset()</code> method. Thus, the
 * message tags fall in the range <I>T</I> .. <I>K</I>&minus;1+<I>T</I>, where
 * <I>K</I> is the total number of workers in all the processes. The program
 * should not use message tags in this range except to send and receive the
 * messages described above.
 * <P>
 * Note that each worker team thread actually creates its own instance of the
 * worker for loop class and passes that instance to the worker region's
 * <code>execute()</code> method. Thus, any fields declared in the worker for loop
 * class will <I>not</I> be shared by all the workers, but instead will be
 * private to each worker.
 * <P>
 * The <code>start()</code> method is intended for performing per-thread
 * initialization before starting the loop iterations. If no such initialization
 * is needed, omit the <code>start()</code> method.
 * <P>
 * The <code>run()</code> method contains the code for the loop. The first and last
 * indexes for a chunk of loop iterations are passed in as arguments. The loop
 * stride, which is always positive, is also explicitly specified as an
 * argument. The worker for loop's <code>run()</code> method must be coded this way:
 * <PRE>
 *     public void run (long first, long last, long stride)
 *         {
 *         for (long i = first; i &lt;= last; i += stride)
 *             {
 *             // Loop body code
 *             . . .
 *             }
 *         }
 * </PRE> with the loop indexes running from <code>first</code> to <code>last</code>
 * inclusive and increasing by <code>stride</code> on each iteration.
 * <P>
 * The <code>finish()</code> method is intended for performing per-thread
 * finalization after finishing the loop iterations. If no such finalization is
 * needed, omit the <code>finish()</code> method.
 * <P>
 * If the worker for loop's <code>start()</code>, <code>run()</code>, or
 * <code>finish()</code> method throws an exception in one of the worker threads,
 * then that worker thread executes no further code in the loop, and the worker
 * region's <code>execute()</code> method throws that same exception in that thread.
 * However, the other worker threads in the worker team continue to execute.
 *
 * @author Alan Kaminsky
 * @version 27-Jan-2010
 */
public abstract class WorkerLongStrideForLoop
        extends WorkerForLoop {

// Exported constructors.
    /**
     * Construct a new worker for loop.
     */
    public WorkerLongStrideForLoop() {
        super();
    }

// Exported operations.
    /**
     * Determine this worker for loop's schedule. Called by the master and
     * worker threads. The schedule determines how the loop iterations are
     * apportioned among the worker team threads. For further information, see
     * class {@linkplain LongSchedule}.
     * <P>
     * The <code>schedule()</code> method may be overridden in a subclass to return
     * the desired schedule. If not overridden, the default is a runtime
     * schedule (see {@link edu.rit.pj.LongSchedule#runtime()}).
     *
     * @return Schedule for this worker for loop.
     */
    public LongSchedule schedule() {
        return LongSchedule.runtime();
    }

    /**
     * Perform per-thread initialization actions before starting the loop
     * iterations. Called by a worker thread.
     * <P>
     * The <code>start()</code> method may be overridden in a subclass. If not
     * overridden, the <code>start()</code> method does nothing.
     *
     * @exception Exception The <code>start()</code> method may throw any exception.
     * @throws java.lang.Exception if any.
     */
    public void start()
            throws Exception {
    }

    /**
     * Send additional input data associated with a task. Called by the master
     * thread. The task is denoted by the given chunk of loop iterations. The
     * input data must be sent using the given communicator, to the given worker
     * process rank, with the given message tag.
     * <P>
     * The <code>sendTaskInput()</code> method may be overridden in a subclass. If
     * not overridden, the <code>sendTaskInput()</code> method does nothing.
     *
     * @param range Chunk of loop iterations.
     * @param comm Communicator.
     * @param wRank Worker process rank.
     * @param tag Message tag.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void sendTaskInput(LongRange range,
            Comm comm,
            int wRank,
            int tag)
            throws IOException {
    }

    /**
     * Receive additional input data associated with a task. Called by a worker
     * thread. The task is denoted by the given chunk of loop iterations. The
     * input data must be received using the given communicator, from the given
     * master process rank, with the given message tag.
     * <P>
     * The <code>receiveTaskInput()</code> method may be overridden in a subclass.
     * If not overridden, the <code>receiveTaskInput()</code> method does nothing.
     *
     * @param range Chunk of loop iterations.
     * @param comm Communicator.
     * @param mRank Master process rank.
     * @param tag Message tag.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void receiveTaskInput(LongRange range,
            Comm comm,
            int mRank,
            int tag)
            throws IOException {
    }

    /**
     * Execute one chunk of iterations of this worker for loop. Called by a
     * worker thread. The <code>run()</code> method must perform the loop body for
     * indexes <code>first</code> through <code>last</code> inclusive, increasing the
     * loop index by <code>stride</code> after each iteration.
     * <P>
     * The <code>run()</code> method must be overridden in a subclass.
     *
     * @param first First loop index.
     * @param last Last loop index.
     * @param stride Loop index stride, always positive.
     * @exception Exception The <code>run()</code> method may throw any exception.
     * @throws java.lang.Exception if any.
     */
    public abstract void run(long first,
            long last,
            long stride)
            throws Exception;

    /**
     * Send additional output data associated with a task. Called by a worker
     * thread. The task is denoted by the given chunk of loop iterations. The
     * output data must be sent using the given communicator, to the given
     * master process rank, with the given message tag.
     * <P>
     * The <code>sendTaskOutput()</code> method may be overridden in a subclass. If
     * not overridden, the <code>sendTaskOutput()</code> method does nothing.
     *
     * @param range Chunk of loop iterations.
     * @param comm Communicator.
     * @param mRank Master process rank.
     * @param tag Message tag.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void sendTaskOutput(LongRange range,
            Comm comm,
            int mRank,
            int tag)
            throws IOException {
    }

    /**
     * Receive additional output data associated with a task. Called by the
     * master thread. The task is denoted by the given chunk of loop iterations.
     * The output data must be received using the given communicator, from the
     * given worker process rank, with the given message tag.
     * <P>
     * The <code>receiveTaskOutput()</code> method may be overridden in a subclass.
     * If not overridden, the <code>receiveTaskOutput()</code> method does nothing.
     *
     * @param range Chunk of loop iterations.
     * @param comm Communicator.
     * @param wRank Worker process rank.
     * @param tag Message tag.
     * @exception IOException Thrown if an I/O error occurred.
     * @throws java.io.IOException if any.
     */
    public void receiveTaskOutput(LongRange range,
            Comm comm,
            int wRank,
            int tag)
            throws IOException {
    }

    /**
     * Perform per-thread finalization actions after finishing the loop
     * iterations. Called by a worker thread.
     * <P>
     * The <code>finish()</code> method may be overridden in a subclass. If not
     * overridden, the <code>finish()</code> method does nothing.
     *
     * @exception Exception The <code>finish()</code> method may throw any
     * exception.
     * @throws java.lang.Exception if any.
     */
    public void finish()
            throws Exception {
    }

    /**
     * Returns the tag offset for this worker for loop. Each message between the
     * master and worker threads is sent with a message tag equal to
     * <I>W</I>+<I>T</I>, where <I>W</I> is the worker index and <I>T</I> is the
     * tag offset.
     * <P>
     * The <code>tagOffset()</code> method may be overridden in a subclass. If not
     * overridden, the <code>tagOffset()</code> returns a default tag offset of
     * <code>Integer.MIN_VALUE</code>.
     *
     * @return Tag offset.
     */
    public int tagOffset() {
        return Integer.MIN_VALUE;
    }

// Hidden operations.
    /**
     * Execute this worker for loop in the master thread.
     *
     * @param range Loop index range.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    void masterExecute(LongRange range)
            throws IOException {
        LongSchedule sch = schedule();
        if (sch.isFixedSchedule()) {
            masterExecuteFixed(range, sch);
        } else {
            masterExecuteNonFixed(range, sch);
        }
    }

    /**
     * Execute this worker for loop in the master thread with a fixed schedule.
     *
     * @param range Loop index range.
     * @param sch Schedule.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    void masterExecuteFixed(LongRange range,
            LongSchedule sch)
            throws IOException {
        int count = myTeam.count;
        Comm comm = myTeam.comm;

        // Send additional task input to each worker.
        sch.start(count, range);
        for (int w = 0; w < count; ++w) {
            LongRange chunk = sch.next(w);
            if (chunk != null) {
                sendTaskInput(chunk, comm, myTeam.workerRank(w), tagFor(w));
            }
        }

        // Receive additional task output from each worker.
        sch.start(count, range);
        for (int w = 0; w < count; ++w) {
            LongRange chunk = sch.next(w);
            if (chunk != null) {
                receiveTaskOutput(chunk, comm, myTeam.workerRank(w), tagFor(w));
            }
        }
    }

    /**
     * Execute this worker for loop in the master thread with a non-fixed
     * schedule.
     *
     * @param range Loop index range.
     * @param sch Schedule.
     *
     * @exception IOException Thrown if an I/O error occurred.
     */
    void masterExecuteNonFixed(LongRange range,
            LongSchedule sch)
            throws IOException {
        int count = myTeam.count;
        sch.start(count, range);
        int remaining = count;
        ObjectItemBuf<LongRange> buf = ObjectBuf.buffer();
        Range tagRange = new Range(tagFor(0), tagFor(count - 1));
        Comm comm = myTeam.comm;

        // Send initial task to each worker.
        for (int w = 0; w < count; ++w) {
            LongRange chunk = sch.next(w);
            buf.item = chunk;
            buf.reset();
            int r = myTeam.workerRank(w);
            int tag = tagFor(w);
            comm.send(r, tag, buf);
            if (chunk == null) {
                --remaining;
            } else {
                sendTaskInput(chunk, comm, r, tag);
            }
        }

        // Repeatedly receive a response from a worker and send next task to
        // that worker.
        while (remaining > 0) {
            CommStatus status = comm.receive(null, tagRange, buf);
            LongRange chunk = buf.item;
            int r = status.fromRank;
            int tag = status.tag;
            int w = workerFor(tag);
            receiveTaskOutput(chunk, comm, r, tag);
            chunk = sch.next(w);
            buf.item = chunk;
            buf.reset();
            comm.send(r, tag, buf);
            if (chunk == null) {
                --remaining;
            } else {
                sendTaskInput(chunk, comm, r, tag);
            }
        }
    }

    /**
     * Execute this worker for loop in a worker thread.
     *
     * @param w Worker index.
     * @param range Loop index range.
     *
     * @exception Exception This method may throw any exception.
     */
    void workerExecute(int w,
            LongRange range)
            throws Exception {
        LongSchedule sch = schedule();
        if (sch.isFixedSchedule()) {
            sch.start(myTeam.count, range);
            workerExecuteFixed(sch.next(w), w);
        } else {
            workerExecuteNonFixed(w);
        }
    }

    /**
     * Execute this worker for loop in a worker thread using a fixed schedule.
     *
     * @param range Chunk of loop iterations.
     * @param w Worker index.
     *
     * @exception Exception This method may throw any exception.
     */
    void workerExecuteFixed(LongRange range,
            int w)
            throws Exception {
        start();
        if (range != null) {
            Comm comm = myTeam.comm;
            int r = myTeam.masterRank();
            int tag = tagFor(w);
            receiveTaskInput(range, comm, r, tag);
            run(range.lb(), range.ub(), range.stride());
            sendTaskOutput(range, comm, r, tag);
        }
        finish();
    }

    /**
     * Execute this worker for loop in a worker thread using a non-fixed
     * schedule.
     *
     * @param w Worker index.
     *
     * @exception Exception This method may throw any exception.
     */
    void workerExecuteNonFixed(int w)
            throws Exception {
        Comm comm = myTeam.comm;
        int r = myTeam.masterRank();
        int tag = tagFor(w);
        start();
        ObjectItemBuf<LongRange> buf = ObjectBuf.buffer();
        for (;;) {
            comm.receive(r, tag, buf);
            LongRange range = buf.item;
            if (range == null) {
                break;
            }
            receiveTaskInput(range, comm, r, tag);
            run(range.lb(), range.ub(), range.stride());

            // The next two statements constitute a critical section; other
            // workers in this team must not send messages in between these two
            // messages, or the master can deadlock.
            synchronized (myTeam) {
                comm.send(r, tag, buf);
                sendTaskOutput(range, comm, r, tag);
            }
        }
        finish();
    }

    /**
     * Returns the message tag for the given worker index.
     *
     * @param w Worker index.
     *
     * @return Message tag.
     */
    private int tagFor(int w) {
        return w + tagOffset();
    }

    /**
     * Returns the worker index for the given message tag.
     *
     * @param tag Message tag.
     *
     * @return Worker index.
     */
    private int workerFor(int tag) {
        return tag - tagOffset();
    }

}
