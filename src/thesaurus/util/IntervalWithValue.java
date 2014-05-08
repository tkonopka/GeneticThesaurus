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

/**
 * A container class for two integers (start/end positions of an interval) and a
 * String value associated with the interval.
 *
 *
 * @author tkonopka
 */
public class IntervalWithValue {

    private final int start, end;
    private final String value;

    public IntervalWithValue(int start, int end, String value) {
        // when declaring the interval, always put the lower value in start
        this.start = Math.min(start, end);
        this.end = Math.max(start, end);
        this.value = value;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public String getValue() {
        return value;
    }

    @Override
    public String toString() {
        return start + "\t" + end + "\t" + value;
    }
}
