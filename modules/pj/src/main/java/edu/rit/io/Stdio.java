//******************************************************************************
//
// File:    Stdio.java
// Package: edu.rit.io
// Unit:    Class edu.rit.io.Stdio
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
package edu.rit.io;

import java.io.InputStream;
import java.io.PrintStream;

/**
 * Class Stdio provides standard I/O streams that can be redirected on a
 * per-thread basis.
 *
 * @author Alan Kaminsky
 * @version 08-Oct-2010
 */
public class Stdio {

// Prevent construction.
    private Stdio() {
    }

// Hidden data members.
    private static ThreadLocal<InputStream> in
            = new ThreadLocal<InputStream>() {
                protected InputStream initialValue() {
                    return System.in;
                }
            };

    private static ThreadLocal<PrintStream> out
            = new ThreadLocal<PrintStream>() {
                protected PrintStream initialValue() {
                    return System.out;
                }
            };

    private static ThreadLocal<PrintStream> err
            = new ThreadLocal<PrintStream>() {
                protected PrintStream initialValue() {
                    return System.err;
                }
            };

// Exported operations.
    /**
     * Get the standard input stream for the calling thread.
     *
     * @return Standard input stream.
     */
    public static InputStream in() {
        return in.get();
    }

    /**
     * Get the standard output stream for the calling thread.
     *
     * @return Standard output stream.
     */
    public static PrintStream out() {
        return out.get();
    }

    /**
     * Get the standard error stream for the calling thread.
     *
     * @return Standard error stream.
     */
    public static PrintStream err() {
        return err.get();
    }

    /**
     * Set the standard input stream for the calling thread. If not set, the
     * default is <code>System.in</code>.
     *
     * @param stream Standard input stream.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>stream</code> is null.
     */
    public static void in(InputStream stream) {
        if (stream == null) {
            throw new NullPointerException("Stdio.in(): stream is null");
        }
        in.set(stream);
    }

    /**
     * Set the standard output stream for the calling thread. If not set, the
     * default is <code>System.out</code>.
     *
     * @param stream Standard output stream.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>stream</code> is null.
     */
    public static void out(PrintStream stream) {
        if (stream == null) {
            throw new NullPointerException("Stdio.out(): stream is null");
        }
        out.set(stream);
    }

    /**
     * Set the standard error stream for the calling thread. If not set, the
     * default is <code>System.err</code>.
     *
     * @param stream Standard error stream.
     * @exception NullPointerException (unchecked exception) Thrown if
     * <code>stream</code> is null.
     */
    public static void err(PrintStream stream) {
        if (stream == null) {
            throw new NullPointerException("Stdio.err(): stream is null");
        }
        err.set(stream);
    }

}
