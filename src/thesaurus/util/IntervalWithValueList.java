/*
 * Copyright 2013 Tomasz Konopka.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package thesaurus.util;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

/**
 *
 * Although it is called a list. this class does not actually implement the List
 * interface. i.e. only a select subset of List-related functionality is
 * implemented here.
 *
 * The list stores multiple intervals and can be sorted. After sorting, the
 * getValue function can be invoked to find what named intervals contain a given
 * position.
 *
 * WARNING: getValue without prior sorting will give rubbish results.
 *
 *
 * @author tkonopka
 */
public class IntervalWithValueList {

    private final ArrayList<IntervalWithValue> list = new ArrayList<IntervalWithValue>();
    private int[] maxend = null;
    private final IntervalWithValueStartComparator startcomparator = new IntervalWithValueStartComparator();

    /**
     * A comparator that will sort the intervals by start position.
     */
    class IntervalWithValueStartComparator implements Comparator {

        @Override
        public int compare(Object o1, Object o2) {
            IntervalWithValue interval1 = (IntervalWithValue) o1;
            IntervalWithValue interval2 = (IntervalWithValue) o2;
            if (interval1.getStart() < interval2.getStart()) {
                return -1;
            } else if (interval1.getStart() > interval2.getStart()) {
                return 1;
            } else {
                return 0;
            }
        }
    }

    /**
     * A comparator that will sort the intervals by start position.
     */
    class IntervalWithValueEndComparator implements Comparator {

        @Override
        public int compare(Object o1, Object o2) {
            IntervalWithValue interval1 = (IntervalWithValue) o1;
            IntervalWithValue interval2 = (IntervalWithValue) o2;
            if (interval1.getEnd() < interval2.getEnd()) {
                return -1;
            } else if (interval1.getEnd() > interval2.getEnd()) {
                return 1;
            } else {
                return 0;
            }
        }
    }

    /**
     * initialize a new list.
     *
     */
    public IntervalWithValueList() {
    }

    /**
     *
     * record an interval into the list.
     *
     * @param start
     *
     * interval start position (included)
     *
     * @param end
     *
     * interval end position (included)
     *
     * @param value
     *
     * some value to associate to the interval
     *
     */
    public void add(int start, int end, String value) {
        IntervalWithValue iwv = new IntervalWithValue(start, end, value);
        list.add(iwv);
    }

    /**
     * rearranges the values within the internal list. This MUST be called
     * before getValues(int position) can output some reasonable answers.
     *
     * Invisibly from the user, this also creates a helper array maxend which
     * helps to speed up the getValues process.
     *
     */
    public void sort() {
        // sort the intervals by start position
        Collections.sort(list, startcomparator);

        // create an array indicating the maximal end position 
        // The size of this array will always be one longer than the interval array
        int listsize = list.size();
        if (listsize > 0) {
            maxend = new int[listsize + 1];
            maxend[0] = list.get(0).getEnd();
            for (int i = 1; i < listsize; i++) {
                maxend[i] = Math.max(maxend[i - 1], list.get(i).getEnd());
            }
            maxend[listsize] = list.get(listsize - 1).getEnd();
        } else {
            maxend = null;
        }
    }

    /**
     *
     * WARNING: the list must have been sorted via sort() for this work
     * properly.
     *
     * @param position
     * @return
     */
    public synchronized ArrayList<String> getValues(int position) {
        //System.out.println("getting for "+position);
        ArrayList<String> ans = new ArrayList<String>();

        if (maxend == null || maxend.length != list.size() + 1) {
            return null;
        }

        // create a dummy interval for this position. This is used for the binary search below.
        IntervalWithValue posinterval = new IntervalWithValue(position, position, null);

        // find out where in the list to start the iteration
        // this is achieved by binary search with the custom comparator, which will
        // only look at the start position of the interval.                        
        int nowindex = Collections.binarySearch(list, posinterval, startcomparator);
        if (nowindex < 0) {
            nowindex = -nowindex - 1;
        }

        // backtrack the index to capture all relevant intervals
        while (nowindex > 0 && maxend[nowindex] >= position) {
            nowindex--;
        }

        while (nowindex < list.size()) {
            IntervalWithValue nowinterval = list.get(nowindex);
            // make sure the position is within the interval
            if (position < nowinterval.getStart()) {
                break;
            }
            // check that the ends of the interval covers the position
            if (position <= nowinterval.getEnd()) {
                ans.add(nowinterval.getValue());
            }
            nowindex++;
        }

        return ans;
    }

    public int size() {
        return list.size();
    }

    public IntervalWithValue get(int index) {
        return list.get(index);
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < size(); i++) {
            IntervalWithValue iwv = list.get(i);
            sb.append(iwv.toString()).append("\n");
        }

        return sb.toString();
    }
}
