This short vignette demonstrates how to use the RGGV package to create
both interactive cluster maps of SNP allele frequencies and
environmental choropleth maps of the type seen in ...

Vignette Info
-------------

Note the various macros within the `vignette` section of the metadata
block above. These are required in order to instruct R how to build the
vignette. Note that you should change the `title` field and the
`\VignetteIndexEntry` to match the title of your vignette.

Styles
------

The `html_vignette` template includes a basic CSS theme. To override
this theme you can specify your own CSS in the document metadata as
follows:

Figures
-------

A GGV interactive cluster map can be generated by calling the ggv
function on a snp rs ID or physical position.

    library(RGGV)

    ggv()

Try zooming in and out on the map, hovering over the pie charts, and
clicking the pie charts.

You can enable figure captions by `fig_caption: yes` in YAML:

    output:
      rmarkdown::html_vignette:
        fig_caption: yes

Then you can use the chunk option `fig.cap = "Your figure caption."` in
**knitr**.

More Examples
-------------

You can write math expressions, e.g. *Y* = *X**β* + *ϵ*, footnotes[1],
and tables, e.g. using `knitr::kable()`.

<table>
<thead>
<tr class="header">
<th></th>
<th align="right">mpg</th>
<th align="right">cyl</th>
<th align="right">disp</th>
<th align="right">hp</th>
<th align="right">drat</th>
<th align="right">wt</th>
<th align="right">qsec</th>
<th align="right">vs</th>
<th align="right">am</th>
<th align="right">gear</th>
<th align="right">carb</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Mazda RX4</td>
<td align="right">21.0</td>
<td align="right">6</td>
<td align="right">160.0</td>
<td align="right">110</td>
<td align="right">3.90</td>
<td align="right">2.620</td>
<td align="right">16.46</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">4</td>
<td align="right">4</td>
</tr>
<tr class="even">
<td>Mazda RX4 Wag</td>
<td align="right">21.0</td>
<td align="right">6</td>
<td align="right">160.0</td>
<td align="right">110</td>
<td align="right">3.90</td>
<td align="right">2.875</td>
<td align="right">17.02</td>
<td align="right">0</td>
<td align="right">1</td>
<td align="right">4</td>
<td align="right">4</td>
</tr>
<tr class="odd">
<td>Datsun 710</td>
<td align="right">22.8</td>
<td align="right">4</td>
<td align="right">108.0</td>
<td align="right">93</td>
<td align="right">3.85</td>
<td align="right">2.320</td>
<td align="right">18.61</td>
<td align="right">1</td>
<td align="right">1</td>
<td align="right">4</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td>Hornet 4 Drive</td>
<td align="right">21.4</td>
<td align="right">6</td>
<td align="right">258.0</td>
<td align="right">110</td>
<td align="right">3.08</td>
<td align="right">3.215</td>
<td align="right">19.44</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">1</td>
</tr>
<tr class="odd">
<td>Hornet Sportabout</td>
<td align="right">18.7</td>
<td align="right">8</td>
<td align="right">360.0</td>
<td align="right">175</td>
<td align="right">3.15</td>
<td align="right">3.440</td>
<td align="right">17.02</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">2</td>
</tr>
<tr class="even">
<td>Valiant</td>
<td align="right">18.1</td>
<td align="right">6</td>
<td align="right">225.0</td>
<td align="right">105</td>
<td align="right">2.76</td>
<td align="right">3.460</td>
<td align="right">20.22</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">1</td>
</tr>
<tr class="odd">
<td>Duster 360</td>
<td align="right">14.3</td>
<td align="right">8</td>
<td align="right">360.0</td>
<td align="right">245</td>
<td align="right">3.21</td>
<td align="right">3.570</td>
<td align="right">15.84</td>
<td align="right">0</td>
<td align="right">0</td>
<td align="right">3</td>
<td align="right">4</td>
</tr>
<tr class="even">
<td>Merc 240D</td>
<td align="right">24.4</td>
<td align="right">4</td>
<td align="right">146.7</td>
<td align="right">62</td>
<td align="right">3.69</td>
<td align="right">3.190</td>
<td align="right">20.00</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">4</td>
<td align="right">2</td>
</tr>
<tr class="odd">
<td>Merc 230</td>
<td align="right">22.8</td>
<td align="right">4</td>
<td align="right">140.8</td>
<td align="right">95</td>
<td align="right">3.92</td>
<td align="right">3.150</td>
<td align="right">22.90</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">4</td>
<td align="right">2</td>
</tr>
<tr class="even">
<td>Merc 280</td>
<td align="right">19.2</td>
<td align="right">6</td>
<td align="right">167.6</td>
<td align="right">123</td>
<td align="right">3.92</td>
<td align="right">3.440</td>
<td align="right">18.30</td>
<td align="right">1</td>
<td align="right">0</td>
<td align="right">4</td>
<td align="right">4</td>
</tr>
</tbody>
</table>

Also a quote using `>`:

> "He who gives up \[code\] safety for \[code\] speed deserves neither."
> ([via](https://twitter.com/hadleywickham/status/504368538874703872))

[1] A footnote here.