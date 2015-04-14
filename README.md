GeneticThesaurus
================

GeneticThesaurus is designed for analysis of genetic variation in repetitive regions. It contains programs that build a table describing genomic regions with similar sequence. It also provides tools that use the resource to annotate variants.

To run the GeneticThesaurus, you will need java (version 6 or above). Execute the program using the command

	java -jar GeneticThesaurus.jar

You may need to provide the java virtual machine access to a larger pool of memory than is allowed by default. In such situations, additionally use the -Xmx flag. For example, to provide 24GB of heap memory, use

	java -Xmx24g -jar GeneticThesaurus.jar



Development and Distribution
----------------------------

The project is developed on github.com and distributed on sourceforge.net under the name "GeneticThesaurus".



Documentation
-------------

Usage information is available on the wiki pages at www.sourceforge.net.



Authors
-------

GeneticThesaurus is described in an academic publication. If you find  GeneticThesaurus useful, please cite

C. Kerzendorfer, T. Konopka, S. Nijman. "A thesaurus of genetic variation for interrogation of repetitive genomic regions" Nucleic Acids Research

T.Konopka, S. Nijman. "GeneticThesaurus: Comparison of matched samples with thesaurus annotation" (under review)

Source code is developed by Tomasz Konopka.



Source code licence
-------------------

Copyright 2013-2015 Tomasz Konopka

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

