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
<link href="https://gitcdn.github.io/bootstrap-toggle/2.2.2/css/bootstrap-toggle.min.css" rel="stylesheet">
<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
<script src="https://code.highcharts.com/highcharts.js"></script>
<script src="https://code.highcharts.com/modules/exporting.js"></script>
<script src="https://code.highcharts.com/modules/export-data.js"></script>
<script src="https://code.highcharts.com/modules/drilldown.js"></script>
<script src="https://gitcdn.github.io/bootstrap-toggle/2.2.2/js/bootstrap-toggle.min.js"></script>
<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;border-color:#aabcfe;}
.tg td{font-family:Arial, sans-serif;font-size:14px;padding:5px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#aabcfe;color:#669;background-color:#e8edff;}
.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#aabcfe;color:#039;background-color:#b9c9fe;}
.tg .tg-hmp3 {
    background-color:#D2E4FC;
    text-align:center;
    vertical-align:bottom;
}
.tg .tg-baqh{text-align:center;vertical-align:top}
.tg .tg-mb3i{background-color:#D2E4FC;text-align:right;vertical-align:top}
.tg .tg-lqy6{text-align:right;vertical-align:top}
.tg .tg-0lax{text-align:left;vertical-align:top}

	
.tg .tg-phtq{background-color:#D2E4FC;border-color:inherit;text-align:left;vertical-align:top}
.tg .tg-c3ow{border-color:inherit;text-align:center;vertical-align:top}
.tg .tg-0pky{border-color:inherit;text-align:left;vertical-align:top}


.border-right {
border-right: solid 1px #039 !important;
}

.tg-hmp3 a {
    cursor: pointer;
}
.count-table .tg-hmp3{
	width: 150px;
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
width: 75%;
margin-left: 12.5%;
/*padding-top: 10px;
*/background-color: #ffffff;
}
.content-body {
padding: 10px;
}
#header {
	text-align: center;
  /*margin-top: 10px;*/
	padding-top: 20px;
	background-color: #f1f1f1;
  font-size: 20px;
}


.read_section  {
	margin-left: 10px
	margin-bottom: 25px;

}

.divtable {
	margin-left: 10px;
	margin-top: 10px;
    overflow : auto;
}

body {
	background-color: #03233f;
    font-size: 16px;
}

.tg .tg-cap {
        border-color:inherit;
        text-align:center;
        vertical-align:top
}

caption {
    display: table-caption;
    text-align: left;
}
/* The navigation bar */
#navbard {
  overflow: hidden;
  position: fixed; /* Set the navbar to fixed position */
  top: 0; /* Position the navbar at the top of the page */
  width: 100%;
  font-size: 15px;

}
</style>
''')
    
## main content div template for each case
contentTemplate = Template(
'''
{% for readLen,count_pLays in counterRes.items()  %}
  {% for lay,counters in count_pLays.items()  %}
    <div id='read_{{readLen}}_{{lay}}' class='read_section'>
    <span>Mapping summary for read length {{readLen}} and {{lay}} layout</span>
	<hr class="hrleft">
    
	<div class="summary-table divtable">
		<table class="tg">
            <caption>Overall mapping statistics</caption>
          <!--  <tr>
    			<th class="tg-cap" colspan="3">Overall mapping statistics</th>
  			</tr>
              -->
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
		    <th scope='row' class="tg-0lax">Total Reads<br></th>
		    <td class="tg-0lax">{{counters.getTotalReads()}}</td>
		    <td class="tg-0lax">--</td>
		  </tr>
		</table>
	</div>
    <br>
    <div class="count-table divtable">
    
		<table class="tg">
        <caption>Correct and Incorrect mapping statistics</caption>
		  <tr>
		    <th class="tg-c3ow" rowspan="3">Genome</th>
		    <th class="tg-c3ow" colspan="2" > <a data-toggle="tooltip" title="Reads mapped to the correct genome regardless of contigs">Correct</a></th>
		    <th class="tg-c3ow" colspan="3"><a data-toggle="tooltip" title="Reads mapped to wrong genomes">InCorrect</a></th>
		  </tr>
		  <tr>
		    
		    <td class="tg-hmp3"  rowspan="2" ><a  data-toggle="tooltip"  title="Correct Uniquely Mapped Reads">Unique</a></td>
		    <td class="tg-hmp3 border-right"  rowspan="2" ><a  data-toggle="tooltip"  title="Correct Multi Mapped Reads">Multi</a></td>
		    <!-- <td class="tg-hmp3" rowspan="2" ><a  data-toggle="tooltip"  title="Incorrect Unique Mapped Reads. Reads mapped to the wrong contig">Unique</a></td> -->
		    <td class="tg-hmp3" rowspan="2"><a  data-toggle="tooltip"  title="Reads originated for source species, but uniquely mapped to other species">Unique Cross</a></td>
            <td class="tg-hmp3" colspan="2"><a  data-toggle="tooltip"  title="Incorrect Multi Mapped Reads" >Multi</a></td>

		  </tr>
          <tr>
           <td class="tg-hmp3"><a  data-toggle="tooltip"  title="Multi Mapped Reads mapped to the source species and other species">Source and Other</a></td>
		    <td class="tg-hmp3" ><a  data-toggle="tooltip"  title="Multi Mapped Reads that did not map to the source species but mapped to other species">Only Other</a></td></tr>
          {% for sp,spId in counters.speciesIds.items() %}
           {% set rowStyle = 'tg-baqh-odd' %} 
           {% set rowHeadStyle = 'tg-0pky' %}
           {% if loop.index is divisibleby 2 %}
               {% set rowStyle= 'tg-baqh-even' %}
               {% set rowHeadStyle = 'tg-phtq' %}
           {% endif %}
           <tr>
		     <th scope='row' class="{{rowHeadStyle}}">{{sp}}<br></th>
             {%if showPercent%}
             <td class="{{rowStyle}}">{{ ( 100*((counters.unique[spId][1]+counters.unique[spId][0])/counters.getTotalReads()))|round(2)}}%</td>
             <td class="border-right {{rowStyle}}">{{( 100*(counters.multiReads[spId][0]/counters.getTotalReads()))|round(2)}}%</td>  
             <!-- <td class="{{rowStyle}}">{{( 100*( counters.unique[spId][1]/counters.getTotalReads()))|round(2)}}%</td> -->
             <td class="{{rowStyle}}">{{( 100*( counters.unique[spId][2]/counters.getTotalReads() ))|round(2)}}%</td>
             <td class="{{rowStyle}}">{{( 100*( counters.multiReads[spId][1]/counters.getTotalReads()))|round(2)}}%</td> 
             <td class="{{rowStyle}}">{{( 100*( counters.multiReads[spId][2]/counters.getTotalReads()))|round(2)}}%</td>
             {%else%}
             <td class="{{rowStyle}}">{{(counters.unique[spId][0]+counters.unique[spId][1])}}</td>
             <td class="border-right {{rowStyle}}">{{counters.multiReads[spId][0]}}</td>  
             <!-- <td class="{{rowStyle}}">{{counters.unique[spId][1]}}</td> -->
             <td class="{{rowStyle}}">{{counters.unique[spId][2]}}</td>
             <td class="{{rowStyle}}">{{counters.multiReads[spId][1]}}</td> 
             <td class="{{rowStyle}}">{{counters.multiReads[spId][2]}}</td>
             {%endif%}
            </tr>
           {% endfor %}
         </table>
    
    </div>
    
    <br>
    <div class='cross-table divtable'>
    {%set nSp = (counters.speciesIds|length) %}
    	<table class="tg">
        <caption>Breakdown of Crossmapping by organisms</caption>

		  <tr>
		    <th class="tg-c3ow" rowspan="2"></th>
		    <th class="tg-c3ow" colspan="{{counters.speciesIds|length}}">Unique Cross</th>
		    <th class="tg-c3ow" colspan="{{counters.speciesIds|length}}">Multi Cross</th>
		    <th class="tg-c3ow" colspan="{{counters.speciesIds|length}}">Total Cross</th>
		  </tr>
          
          <tr>
		 
		    {% for repI in range(0,3) %}
                {% for sp in counters.speciesIds %}
                    {% set clmStyle= 'clm' %}
                    {% if loop.index ==  (counters.speciesIds|length) %}
                         {% set clmStyle= 'border-right' %}
                    {% endif %}
                    <td class="tg-hmp3 {{clmStyle}}">{{sp}}</td>
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
                     {% set clmStyle = 'clm' %}
                     {% if loop.index ==  (counters.speciesIds|length) %}
                         {% set clmStyle= 'border-right' %}
                     {% endif %}
                  <td class="{{rowStyle}} {{clmStyle}}">
                  {%if sp==tar_sp%}
                  -
                  {%else%}
                  <a data-toggle="tooltip" data-html="true"  title="{% if showPercent  %}{{counters.getPerSpCrossMapped(sp,tar_sp,'Unique')}} #Total Mapped Reads{%else%}{{counters.getPerSpCrossMapped(sp,tar_sp,'Unique',percent = True)}}% of Total Mapped Reads{%endif%}.<br> {{counters.getPerSpCrossMapped(sp,tar_sp,'Unique',percent = True , relative = True)}}% of Total Spp. Mapped Reads. ">
                  {{counters.getPerSpCrossMapped(sp,tar_sp,'Unique' , percent = showPercent )}}{% if showPercent  %}%{%endif%}
                  </a>
                  {%endif%}
                  </td>
                 {%endfor%}
                 {% for tar_sp in counters.speciesIds %}
                     {% set clmStyle= 'clm' %}
                     {% if loop.index ==  (counters.speciesIds|length) %}
                         {% set clmStyle= 'border-right' %}
                     {% endif %}
                  <td class="{{rowStyle}} {{clmStyle}}">
                  {%if sp==tar_sp%}
                  -
                  {%else%}
                  <a data-toggle="tooltip" data-html="true"  title="{% if showPercent  %}{{counters.getPerSpCrossMapped(sp,tar_sp,'Multi')}} #Total Mapped Reads{%else%}{{counters.getPerSpCrossMapped(sp,tar_sp,'Multi',percent = True)}}% of Total Mapped Reads{%endif%}.<br> {{counters.getPerSpCrossMapped(sp,tar_sp,'Multi',percent = True , relative = True)}}% of Total Spp. Mapped Reads. ">
                  {{counters.getPerSpCrossMapped(sp,tar_sp,'Multi' , percent = showPercent )}}{% if showPercent  %}%{%endif%}
                  </a>
                  {%endif%}
                  </td>
                 {%endfor%}
                 {% for tar_sp in counters.speciesIds %}
                     {% set clmStyle= 'clm' %}
                     {% if loop.index == (counters.speciesIds|length) %}
                         {% set clmStyle= 'border-right' %}
                     {% endif %}
                  <td class="{{rowStyle}} {{clmStyle}}">
                  {%if sp==tar_sp%}
                  -
                  {%else%}
                       <a data-toggle="tooltip" data-html="true"  title="{% if showPercent  %}{{counters.getPerSpCrossMapped(sp,tar_sp)}} #Total Mapped Reads{%else%}{{counters.getPerSpCrossMapped(sp,tar_sp,percent = True)}}% of Total Mapped Reads{%endif%}.<br> {{counters.getPerSpCrossMapped(sp,tar_sp,percent = True , relative = True)}}% of Total Spp. Mapped Reads. ">
                       {{counters.getPerSpCrossMapped(sp,tar_sp , percent = showPercent )}}{% if showPercent  %}%{%endif%}
                       </a>
                 {%endif%}
                 </td>
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
        text: 'Details of crossmapped reads per species'
    },
    subtitle: {
        text: 'Click the columns to view details of cross mapped count. Click read lengths at the bottom to add/remove barplots'
    },
    xAxis: {
        type: 'category'
    },
    yAxis: {
        title: {
            text: '{% if showPercent  %}% of {%else%}#{%endif%} Wrong Mapped Reads'
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
                format: '{point.y:.1f}{% if showPercent %}%{%endif%}'
            }
        }
    },

    tooltip: {
        headerFormat: '<span style="font-size:11px">{series.name}</span><br>',
        pointFormat: '<span style="color:{point.color}">{point.name}</span>: <b>{point.y:.3f}{% if showPercent %}%{%endif%}    </b> Total Cross-Mapped Reads<br/>'
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
                            "y" : {{counters.getTotalCrossMapped(readLen,lay,sp , percent = showPercent)}},
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
    
    
### TODO :: try to replace "name" : "{{sp}}"  with "name" : {{readLen}}_{{lay}}
drilldownTemplate = Template(
'''
"series" : [
        {% for readLen,count_pLays in counterRes.items()  %}
        {% for lay,counters in count_pLays.items()  %}
            {% for sp,id in counters.speciesIds.items() %}
            {                    
                "name" : "{{sp}} {{readLen}} {{lay}}", 
                "id" : "{{sp}}_{{readLen}}_{{lay}}",
                "data" : 
                    [
                        {%for spVPair in counters.getSpeciesCorssMapped(sp, percent = showPercent)%}
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

 //   subtitle: {
 //       text: '...'
 //   },
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
            text: '{% if showPercent  %}% of {%else%}#{%endif%}Total reads'
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

    series: [
    {%if not args.read_layout == 'PE' %}
    {
        name: 'Total SE',
        data: 
            [
            {% for readLen,count_pLays in counterRes.items()  %}
                {{count_pLays['SE'].getReadLenSeries('Total',percent = showPercent)}},
            {% endfor %}
            ]
    }, 
    {
        name: 'Unique SE',
        data:
            [
            {% for readLen,count_pLays in counterRes.items()  %}
                {{count_pLays['SE'].getReadLenSeries('Unique',percent = showPercent)}},
            {% endfor %}
            ]
    }, 
    {
        name: 'Multi SE',
        data: 
            [
            {% for readLen,count_pLays in counterRes.items()  %}
                {{count_pLays['SE'].getReadLenSeries('Multi',percent = showPercent)}},
            
            {% endfor %}
            ]
    },
    {%endif%}
    {%if not args.read_layout == 'SE' %}
    {
        name: 'Total PE',
       data: 
           [
            {% for readLen,count_pLays in counterRes.items()  %}
                {{count_pLays['PE'].getReadLenSeries('Total',percent = showPercent)}},
            {% endfor %}
            ]
    }, {
        name: 'Unique PE',
        data:
            [
            {% for readLen,count_pLays in counterRes.items()  %}
                {{count_pLays['PE'].getReadLenSeries('Unique',percent = showPercent)}},
            {% endfor %}
            ]
    },{
        name: 'Multi PE',
        data:
            [
            {% for readLen,count_pLays in counterRes.items()  %}
            
                {{count_pLays['PE'].getReadLenSeries('Multi',percent = showPercent)}},
            {% endfor %}
            ]
    }
        {%endif%}
    ],

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

barchart2DivtableTemplate = Template('''
<table style="table-layout: fixed;width: 100%">
  


{% for readLen,count_pLays in counterRes.items()  %}
<tr>
    {% for lay,counters in count_pLays.items()  %}
    <td> <div id="{{lay}}_{{readLen}}" class="divchart" style="width:100%"> </div> </td>
    {%endfor%}
<tr>
{%endfor%}
</table>                     
                                    ''')
   
barchart2Template = Template('''
  
{% for readLen,count_pLays in counterRes.items()  %}
    {% for lay,counters in count_pLays.items()  %}
Highcharts.chart('{{lay}}_{{readLen}}', {
    chart: {
        type: 'column'
    },
    title: {
        text: 'Total Cross-Mapped Reads for Read length {{readLen}} {{lay}}'
    },
    xAxis: {
        categories: [ {% for sp,id in counters.speciesIds.items() %} '{{sp}}',  {%endfor%} ]
    },
    yAxis: {
        min: 0,
        title: {
            text: '{% if showPercent  %}% of {%else%}#{%endif%}Wrong Mapped Reads'
        },
        stackLabels: {
            enabled: true,
            style: {
                fontWeight: 'bold',
                color: (Highcharts.theme && Highcharts.theme.textColor) || 'gray'
            }
        }
    },
    legend: {
        align: 'right',
        x: -30,
        verticalAlign: 'top',
        y: 25,
        floating: true,
        backgroundColor: (Highcharts.theme && Highcharts.theme.background2) || 'white',
        borderColor: '#CCC',
        borderWidth: 1,
        shadow: false
    },
    tooltip: {
        headerFormat: '<b>{point.x}</b><br/>',
        pointFormat: '{series.name}: {point.y}{% if showPercent  %}%{%endif%}<br/>Total: {point.stackTotal}{% if showPercent  %}%{%endif%}'
    },
    plotOptions: {
        column: {
            stacking: 'normal',
            dataLabels: {
                enabled: true,
                color: (Highcharts.theme && Highcharts.theme.dataLabelsColor) || 'white'
            }
        }
    },
    series: [
    {% for sp,id in counters.speciesIds.items() %} 
    { 
     name: '{{sp}}',
     data : [
         {% for tar_sp in counters.speciesIds %}
         {{counters.getPerSpCrossMapped(tar_sp,sp,percent = showPercent)}} ,
         {%endfor%}
         ]
    },
    {%endfor%} 
    ]
    
    
});
    
    
    
    
    {%endfor%}

{%endfor%}                     
                                    ''')

 
reportTemplete = Template('''
<!doctype html>
<html>
<head>
    {% include headTemplate %}
</head>
<body>
<div id="navbard">
    <div style="float: right;color: white;"> Switch to 
    <input id="toggle-percent" type="checkbox" data-toggle="toggle" data-on="Count" data-off="Percent " data-onstyle="warning" data-offstyle="info" size='small'>
    </div>
</div>
<div id='content'>
	<div id='header'>
		Crossmapper Summary Report
		<hr style="width: 100%">
	</div>
    
    
    
    <div class="total-summary " style="width:100%;">
		<div id="container1" style="width:90%; margin-left:5%;"></div>

		<div id="barchartcontainer" style="width:90%; margin-left:5%;">
        {% if not args.groupBarChart %}
        {% include barchart2DivtableTemplate %}
        {%endif%}
		</div>

   </div>
   <hr class='hrsep'>
    
    <div   id="content-count" class='content-body'>
    {% set showPercent= False %}
    {% include contentTemplate %}
    </div>

    <div   id="content-percent" class='content-body'>
    {% set showPercent= True %}
    {% include contentTemplate %}
    </div>



<div>

<script type="text/javascript">
$(document).ready(function(){
    
    $('[data-toggle="tooltip"]').tooltip();   
    Highcharts.setOptions({
        lang: {
            drillUpText: '<< Back'
        }
    });
     $('#toggle-percent').change(function() {
      
      if ( $(this).prop('checked')) {
        $('#content-count').hide();
        $('#content-percent').show();
        
         {% set showPercent= True %}
        {% include lineChartTemplate %}
        {% if args.groupBarChart  %}
        {% include barGroupChartTemplate %}
        {% else %}
            {% include barchart2Template %}
        {% endif %}
        
      }
      else{
        $('#content-count').show();
        $('#content-percent').hide();
        
        {% set showPercent= False %}
        {% include lineChartTemplate %}
        {% if args.groupBarChart  %}
        {% include barGroupChartTemplate %}
        {% else %}
            {% include barchart2Template %}
        {% endif %}
      }
    

    
    
    })
    
    $('#toggle-percent').prop('checked',false).change()


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
            barchart2DivtableTemplate = barchart2DivtableTemplate,
            barchart2Template =barchart2Template ,
            counterRes = resCounters,
            args = args)
    
    with open(reportFilePath, "w+") as fh:
        fh.write(reportHTML)
    
    return




