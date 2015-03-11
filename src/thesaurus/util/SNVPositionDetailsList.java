/*
 * Copyright 2015 Tomasz Konopka.
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

/**
 * A class with an arraylist of SNVPosition objects. Provides some additional search
 * functions for the array.
 * 
 * @author tkonopka
 */
public class SNVPositionDetailsList extends ArrayList<SNVPositionDetails> {

        public SNVPositionDetailsList(int allsize) {
            super(allsize);
        }

        /**
         *
         * @param gpdfind
         * @param vcomp
         * @return
         *
         * an SNVPositionDetails object stored in this list that matches the
         * coordinates given in the input.
         *
         */
        public SNVPositionDetails find(SNVPosition gpdfind, SNVPositionComparator vcomp) {
            int a = Collections.binarySearch(this, gpdfind, vcomp);
            if (a < 0) {
                return null;
            } else {
                return this.get(a);
            }
        }

        public int findIndex(SNVPosition gpdfind, SNVPositionComparator vcomp) {            
            return Collections.binarySearch(this, gpdfind, vcomp);            
        }
                         
}
