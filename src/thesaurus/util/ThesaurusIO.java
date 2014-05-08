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

/**
 * Miscelleaneous static functions for printing output to the console
 *
 * @author tkonopka
 */
public class ThesaurusIO {

    private static String spaceRepeat(int n) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < n; i++) {
            sb.append(" ");
        }
        return (sb.toString());
    }

    public static void printHelpItem(String key, String description) {
        int keylength = key.length();

        // split the description into multiple lines
        String[] dtokens = description.split(" ");
        ArrayList<String> dlines = new ArrayList<String>();
        String nowline = "";
        int nowtoken = 0;
        while (nowtoken < dtokens.length) {
            if (nowline.length() > 50) {
                dlines.add(nowline + "\n");
                nowline = "";
            }
            nowline += dtokens[nowtoken] + " ";
            nowtoken++;
        }
        dlines.add(nowline + "\n");


        StringBuilder sb = new StringBuilder();
        if (key.length() < 30) {
            sb.append("  ").append(key);
            sb.append(spaceRepeat(26 - keylength));
        } else {
            sb.append("  ").append(key).append("\n");
            sb.append(spaceRepeat(28));
        }
        sb.append("- ").append(dlines.get(0));
        for (int i = 1; i < dlines.size(); i++) {
            sb.append(spaceRepeat(30)).append(dlines.get(i));
        }

        System.out.print(sb.toString());
    }
}
