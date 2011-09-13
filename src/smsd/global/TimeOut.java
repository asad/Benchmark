/**
 *
 * Copyright (C) 2006-2011  Syed Asad Rahman <asad@ebi.ac.uk>
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package smsd.global;

import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;

/**
 * Class that manages MCS timeout.
 * @cdk.module smsd
 * @cdk.githash
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
@TestClass("org.openscience.cdk.smsd.global.TimeOutTest")
public class TimeOut {

    private static TimeOut instance = null;
    private double cdkMCSTimeout = -1;
    private double mcsPlusTimeout = -1;
    private double vfTimeout = -1;
    private boolean timeOutFlag = false;

    /**
     * Get Instance of the timeout. This starts the timeout counter.
     * @return Instance
     */
    @TestMethod("testGetInstance")
    public static synchronized TimeOut getInstance() {
        if (instance == null) {
            // it's ok, we can call this constructor
            instance = new TimeOut();
        }
        return instance;
    }

    protected TimeOut() {
    }

    /**
     * set cutoff value for cdkMCSTimeout out eg. -1 for infinite and 0.23 for
     * 23 seconds.
     * @param timeout
     */
    @TestMethod("testSetTimeOut")
    public synchronized void setCDKMCSTimeOut(double timeout) {
        this.cdkMCSTimeout = timeout;
    }

    /**
     * Return cutoff value for cdkMCSTimeout out.
     * @return cdkMCSTimeout out cutoff value
     */
    @TestMethod("testSetTimeOut")
    public synchronized double getCDKMCSTimeOut() {
        return cdkMCSTimeout;
    }

    /**
     * Return true if its a timeout else return false.
     * @return the timeout flag
     */
    @TestMethod("testIsTimeOutFlag")
    public synchronized boolean isTimeOutFlag() {
        return timeOutFlag;
    }

    /**
     * Set true if timeout occures else false
     * @param timeOut the timeout flag to set
     */
    @TestMethod("testSetTimeOutFlag")
    public synchronized void setTimeOutFlag(boolean timeOut) {
        this.timeOutFlag = timeOut;
    }

    /**
     * @return time out for the MCS Plus algorithm
     */
    public synchronized double getMCSPlusTimeout() {
        return mcsPlusTimeout;
    }

    /**
     * @param mcsPlusTimeout time out for mcsPlus
     */
    public synchronized void setMCSPlusTimeout(double mcsPlusTimeout) {
        this.mcsPlusTimeout = mcsPlusTimeout;
    }

    /**
     * @return time out for the VF algorithm
     */
    public synchronized double getVFTimeout() {
        return vfTimeout;
    }

    /**
     * @param vfTimeout time out for the VF algorithm
     */
    public synchronized void setVFTimeout(double vfTimeout) {
        this.vfTimeout = vfTimeout;
    }
}
