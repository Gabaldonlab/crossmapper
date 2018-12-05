"""
Created on Tue Dec  4 09:24:15 2018

@author: ahafez
"""
from jinja2 import Template
import os
from crossmap.helpers import getLogger


#######################################    report Templates
headTemplate = Template('''
<meta charset="utf-8"/>
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
<script src="https://code.highcharts.com/highcharts.js"></script>
<script src="https://code.highcharts.com/modules/exporting.js"></script>
<script src="https://code.highcharts.com/modules/export-data.js"></script>
<script src="https://code.highcharts.com/modules/drilldown.js"></script>

<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;border-color:#aabcfe;}
.tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#aabcfe;color:#669;background-color:#e8edff;}
.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#aabcfe;color:#039;background-color:#b9c9fe;}
.tg .tg-hmp3{background-color:#D2E4FC;text-align:left;vertical-align:top}
.tg .tg-baqh{text-align:center;vertical-align:top}
.tg .tg-mb3i{background-color:#D2E4FC;text-align:right;vertical-align:top}
.tg .tg-lqy6{text-align:right;vertical-align:top}
.tg .tg-0lax{text-align:left;vertical-align:top}

	
.tg .tg-phtq{background-color:#D2E4FC;border-color:inherit;text-align:left;vertical-align:top}
.tg .tg-c3ow{border-color:inherit;text-align:center;vertical-align:top}
.tg .tg-0pky{border-color:inherit;text-align:left;vertical-align:top}




.tg-hmp3 a {
    cursor: pointer;
}
.count-table .tg-hmp3{
	width: 120px;
}
.tg-baqh-even{
background-color: #D2E4FC !important;
}
.tg .tg-svo0{background-color:#D2E4FC;border-color:inherit;text-align:center;vertical-align:top}
.tg .tg-j0tj{background-color:#D2E4FC;text-align:center;vertical-align:top}
.hrleft {
margin-left: 0;
width: 25%;	
}

.hrsep {
border-top: 4px double #03233f;
}

#content {
width: 70%;
margin-left: 15%;
/*padding-top: 10px;
*/background-color: #ffffff;
}
#content-body {
padding: 10px;
}
#header {
	text-align: center;
	padding-top: 10px;
	background-color: #f1f1f1;
    font-size: 20px;
}


.read_section  {
	margin-left: 10px
	margin-bottom: 25px;

}

.divtable {
	margin-left: 10px;
	margin-top: 10px
}

body {
	background-color: #03233f;
    font-size: 16px;
}
</style>
''')
    
## main content div template for each case
contentTemplate = Template(
'''
{% for readLen,count_pLays in counterRes.items()  %}
  {% for lay,counters in count_pLays.items()  %}
    <div id='read_{{readLen}}_{{lay}}' class='read_section'>
    <span> Counters Summary for read length {{readLen}} and ({{lay}}) layout </span>
	<hr class="hrleft">
    
	<div class="summary-table divtable">
		<table class="tg">
		  <tr>
		    <th class="tg-c3ow"></th>
		    <th class="tg-c3ow"><br>#Reads</th>
		    <th class="tg-0pky"><br>Percent</th>
		  </tr>
		  <tr>
		    <th scope='row' class="tg-phtq">Unique Mapped Reads<br></th>
		    <td class="tg-phtq">{{counters.getTotalUniqueMapped()}}</td>
		    <td class="tg-phtq">{{counters.getTotalUniqueMapped(percent=True)}}%</td>
		  </tr>
		  <tr>
		    <th scope='row' class="tg-0pky">Multi Mapped Reads<br></th>
		    <td class="tg-0pky">{{counters.getTotalMultiMapped()}}</td>
		    <td class="tg-0pky">{{counters.getTotalMultiMapped(percent=True)}}%</td>
		  </tr>
		  <tr>
		    <th scope='row' class="tg-phtq">UnMapped Mapped Reads<br></th>
		    <td class="tg-phtq">{{counters.getTotalUnmapped()}}</td>
		    <td class="tg-phtq">{{counters.getTotalUnmapped(percent=True)}}%</td>
		  </tr>
		  <tr>
		    <th scope='row' class="tg-0lax">Totol Reads<br></th>
		    <td class="tg-0lax">{{counters.getTotalReads()}}</td>
		    <td class="tg-0lax">--</td>
		  </tr>
		</table>
	</div>
    <br>
    <div class="count-table divtable">
    
		<table class="tg">
		  <tr>
		    <th class="tg-c3ow" rowspan="2">Genome</th>
		    <th class="tg-c3ow" colspan="2">Correct</th>
		    <th class="tg-c3ow" colspan="4">UnCorrect</th>
		  </tr>
		  <tr>
		    
		    <td class="tg-hmp3"><a  data-toggle="tooltip"  title="Correct Uniquely Mapped Reads">Unique</a></td>
		    <td class="tg-hmp3"><a  data-toggle="tooltip"  title="Correct Multi Mapped Reads">Multi</a></td>
		    <td class="tg-hmp3"><a  data-toggle="tooltip"  title="Uncorrect Unique Mapped Reads. Reads mapped to the wrong contif">Unique</a></td>
		    <td class="tg-hmp3"><a  data-toggle="tooltip"  title="Reads originated from one speices and mapped uniquely to another species">Unique Cross</a></td>
		    <td class="tg-hmp3"><a  data-toggle="tooltip"  title="Multi Mapped Reads mapped to the source species and other species">Multi_org other</a></td>
		    <td class="tg-hmp3"><a  data-toggle="tooltip"  title="Multi Mapped Reads that did not map to the source species but mapped to other species">Multi_onlyOther</a></td>
		  </tr>
          
          {% for sp,spId in counters.speciesIds.items() %}
           {% set rowStyle = 'tg-baqh-odd' %} 
           {% set rowHeadStyle = 'tg-0pky' %}
           {% if loop.index is divisibleby 2 %}
               {% set rowStyle= 'tg-baqh-even' %}
               {% set rowHeadStyle = 'tg-phtq' %}
           {% endif %}
           <tr>
		     <th scope='row' class="{{rowHeadStyle}}">{{sp}}<br></th>
             <td class="{{rowStyle}}">{{counters.unique[spId][0]}}</td>
             <td class="{{rowStyle}}">{{counters.multiReads[spId][0]}}</td>  
             <td class="{{rowStyle}}">{{counters.unique[spId][1]}}</td>
             <td class="{{rowStyle}}">{{counters.unique[spId][2]}}</td>
             <td class="{{rowStyle}}">{{counters.multiReads[spId][1]}}</td> 
             <td class="{{rowStyle}}">{{counters.multiReads[spId][2]}}</td>
            </tr>
           {% endfor %}
         </table>
    
    </div>
    
    <br>
    <div class='cross-table divtable'>
    	<table class="tg">
		  <tr>
		    <th class="tg-c3ow" rowspan="2"></th>
		    <th class="tg-c3ow" colspan="{{counters.speciesIds|length}}">Unique Cross</th>
		    <th class="tg-c3ow" colspan="{{counters.speciesIds|length}}">Multi Cross</th>
		    <th class="tg-c3ow" colspan="{{counters.speciesIds|length}}">Total Cross</th>
		  </tr>
          
          <tr>
		 
		    {% for repI in range(0,3) %}
                {% for sp in counters.speciesIds %}
                    <td class="tg-hmp3">{{sp}}</td>
                {% endfor %}
            {% endfor %}
		  </tr>
           {% for sp in counters.speciesIds %}
           {% set rowStyle = 'tg-baqh-odd' %} 
           {% set rowHeadStyle = 'tg-0pky' %}
           {% if loop.index is divisibleby 2 %}
               {% set rowStyle= 'tg-baqh-even' %}
               {% set rowHeadStyle = 'tg-phtq' %}
           {% endif %}
           <tr>
		     <th scope='row' class="{{rowHeadStyle}}">{{sp}}<br></th>
                 {% for tar_sp in counters.speciesIds %}
                     {% set clmStyle= 'clm' %}
                     {% if loop.index ==  1 %}
                         {% set clmStyle= 'clm-first' %}
                     {% endif %}
                  <td class="{{rowStyle}} {{clmStyle}}">{{counters.getPerSpCrossMapped(sp,tar_sp,'Unique')}}</td>
                 {%endfor%}
                 {% for tar_sp in counters.speciesIds %}
                     {% set clmStyle= 'clm' %}
                     {% if loop.index ==  1 %}
                         {% set clmStyle= 'clm-first' %}
                     {% endif %}
                  <td class="{{rowStyle}} {{clmStyle}}">{{counters.getPerSpCrossMapped(sp,tar_sp,'Multi')}}</td>
                 {%endfor%}
                 {% for tar_sp in counters.speciesIds %}
                     {% set clmStyle= 'clm' %}
                     {% if loop.index == 1 %}
                         {% set clmStyle= 'clm-first' %}
                     {% endif %}
                  <td class="{{rowStyle}} {{clmStyle}}">{{counters.getPerSpCrossMapped(sp,tar_sp)}}</td>
                 {%endfor%}
            </tr>
           {% endfor %}
          
          
         </table>
    </div>
    
    <hr class='hrsep'>
    </div>
 {% endfor %}
{% endfor %}

''')
####################################################################
barGroupChartTemplate = Template(
'''
Highcharts.chart('barchartcontainer', {
    chart: {
        type: 'column'
    },
    title: {
        text: 'Total/Details Cross-Mapped Reads per Species'
    },
    subtitle: {
        text: 'Click the columns to view details cross mapped count'
    },
    xAxis: {
        type: 'category'
    },
    yAxis: {
        title: {
            text: '#Wrong Mapped Reads'
        }

    },
    legend: {
        enabled: true
    },
    plotOptions: {
        series: {
            borderWidth: 0,
            dataLabels: {
                enabled: true,
                format: '{point.y:.1f}'
            }
        }
    },

    tooltip: {
        headerFormat: '<span style="font-size:11px">{series.name}</span><br>',
        pointFormat: '<span style="color:{point.color}">{point.name}</span>: <b>{point.y:.2f}</b> Total Cross-Mapped Reads<br/>'
    },
    {% include seriesTemplate %},
    "drilldown": { {% include drilldownTemplate %}
    }
    });
''')

## this is for the barGroupChartTemplate
seriesTemplate = Template(
'''
"series" : [
     {% for readLen,count_pLays in counterRes.items()  %}
        {% for lay,counters in count_pLays.items()  %}
        {
                "name": "{{readLen}} {{lay}}",
                "colorByPoint": false,
                "data" : [
                        {% for sp,id in counters.speciesIds.items() %}
                        {
                            "name" : "{{sp}}",
                            "y" : {{counters.getTotalCrossMapped(readLen,lay,sp)}},
                            "drilldown" : "{{sp}}_{{readLen}}_{{lay}}"
                        },
                        {% endfor %}
                ]
        },
        {% endfor %}
    {% endfor %}
]
'''
        )
    
drilldownTemplate = Template(
'''
"series" : [
        {% for readLen,count_pLays in counterRes.items()  %}
        {% for lay,counters in count_pLays.items()  %}
            {% for sp,id in counters.speciesIds.items() %}
            {                    
                "name" : "{{sp}}",
                "id" : "{{sp}}_{{readLen}}_{{lay}}",
                "data" : 
                    [
                        {%for spVPair in counters.getSpeciesCorssMapped(sp)%}
                        [
                            "{{spVPair[0]}}",
                            {{spVPair[1]}}
                        ],
                        {% endfor %}
                    ]
            },
            {% endfor %}
            
        {% endfor %}
    {% endfor %}     
]
''')
##########################################################
    
lineChartTemplate = Template(
'''
Highcharts.chart('container1', {

    title: {
        text: 'Total Cross-Mapped Reads per read length'
    },

    subtitle: {
        text: '...'
    },
	 xAxis: {
	 	 title: {
            text: 'Read length'
        },
        categories: [ 
         {% for readLen in counterRes  %}
             '{{readLen}}',
         {% endfor %}
         ]
    },
    yAxis: {
        title: {
            text: 'Total reads'
        }
    },
    legend: {
        layout: 'vertical',
        align: 'right',
        verticalAlign: 'middle'
    },
    tooltip: {
        shared: true,
        crosshairs: true
    },
    plotOptions: {
        series: {
            label: {
                connectorAllowed: false
            }
        }
    },

    series: [{
        name: 'Total SE',
        data: 
            [
            {% for readLen,count_pLays in counterRes.items()  %}
                {{count_pLays['SE'].getReadLenSeries('Total')}},
            {% endfor %}
            ]
    }, {
        name: 'Unique SE',
        data:
            [
            {% for readLen,count_pLays in counterRes.items()  %}
                {{count_pLays['SE'].getReadLenSeries('Unique')}},
            {% endfor %}
            ]
    }, {
        name: 'Multi SE',
        data: 
            [
            {% for readLen,count_pLays in counterRes.items()  %}
                {{count_pLays['SE'].getReadLenSeries('Multi')}},
            
            {% endfor %}
            ]
    }, {
        name: 'Total PE',
       data: 
           [
            {% for readLen,count_pLays in counterRes.items()  %}
                {{count_pLays['PE'].getReadLenSeries('Total')}},
            {% endfor %}
            ]
    }, {
        name: 'Unique PE',
        data:
            [
            {% for readLen,count_pLays in counterRes.items()  %}
                {{count_pLays['PE'].getReadLenSeries('Unique')}},
            {% endfor %}
            ]
    },{
        name: 'Multi PE',
        data:
            [
            {% for readLen,count_pLays in counterRes.items()  %}
            
                {{count_pLays['PE'].getReadLenSeries('Multi')}},
            {% endfor %}
            ]
    }],

    responsive: {
        rules: [{
            condition: {
                maxWidth: 500
            },
            chartOptions: {
                legend: {
                    layout: 'horizontal',
                    align: 'center',
                    verticalAlign: 'bottom'
                }
            }
        }]
    }

});
''')
    
reportTemplete = Template('''
<!doctype html>
<html>
<head>
    {% include headTemplate %}
</head>
<body>

<div id='content'>
	<div id='header'>
		Crossmapper Summary Report
		<hr style="width: 100%">
	</div>
    
    <div  id='content-body'>
    {% include contentTemplate %}
    </div>


    <div class="total-summary " style="width:100%;">
		<div id="container1" style="width:90%; margin-left:5%;"></div>
		</div> 

		<div id="barchartcontainer" style="width:90%; margin-left:5%;"></div>
		</div>

   </div>
<div>

<script type="text/javascript">
$(document).ready(function(){
    
    $('[data-toggle="tooltip"]').tooltip();   
    {% include lineChartTemplate %}
    {% include barGroupChartTemplate %}
});
</script>

</body>
</html>
''')

def createHTMLReport(resCounters,args):
    logger = getLogger()
    reportFilePath = os.path.join(args.out_dir,"report.html")
    logger.info(f"Creating Report File : {reportFilePath}")
    reportHTML = reportTemplete.render(
            headTemplate = headTemplate ,
            contentTemplate = contentTemplate,
            barGroupChartTemplate = barGroupChartTemplate,
            lineChartTemplate = lineChartTemplate,
            seriesTemplate = seriesTemplate ,
            drilldownTemplate = drilldownTemplate ,
            counterRes = resCounters,
            args = args)
    
    with open(reportFilePath, "w") as fh:
        fh.write(reportHTML)
    
    return




