var coefficient = 10;
var genome_id_list = ["gid_ref", "gid_22", "gid_23", "gid_20", "gid_21", "gid_24", "gid_25", "gid_29", "gid_31", "gid_26", "gid_27", "gid_28", "gid_30", "gid_32"];

$(document).ready(function(){
 
// Draw SNPs Table
    $.ajax({
        async:false, dataType:"json",
        url:"js/SNP_table.json",
        success: function(data){
            dr_snp(data);
        }
    });
    
// Draw Phylogeny Tree
    $.ajax({
        async:false, dataType:"json",
        url:"js/strep-tree.json",
        success: function(data){
            drtree(data);
        }
    });
    
// Draw Reference Genome
    $.ajax({
        async:false, dataType:"json",
        url:"js/ref.json",
        success: function(data){
            drref(data);
        }
    });
});

function drref(df){
    var max = 0;
    var maxloc = 0;
    var min = 5000;
    var minloc = 5000;
    var x_shifted = 10;
    var y_level_scaler = 15.5;
    var y_level_contig = 30.5;
    $.each(df, function(key, value){
        if (value.stop >= max){
            max = value.stop;
        };
        if (value.start <= min){
            min = value.start;
        };
        if (value.stoploc >= maxloc){
            maxloc = value.stoploc;
        };
        if (value.startloc <= minloc){
            minloc = value.startloc;
        };
    });
    
    // Create an svg area
    var svg = d3.select("#refmap")
                .append("svg")
                .attr("width", maxloc/coefficient + 50)
                .attr("height", "40");
    
    // Draw scalers
    var scaler = svg.append("line")
                    .attr("x1", minloc/coefficient + x_shifted)
                    .attr("y1", y_level_scaler)
                    .attr("x2", maxloc/coefficient + x_shifted)
                    .attr("y2", y_level_scaler)
                    .attr("stroke", "black");
    
    var hpl = svg.append("text")
                        .attr("x", minloc/coefficient + x_shifted)
                        .attr("y", y_level_scaler - 5)
                        .attr("font-size", "10")
                        .attr("font-family", "sans-serif")
                        .text("Chromosome")
    
    var tpl = svg.append("text")
                        .attr("x", maxloc/coefficient + 2)
                        .attr("y", y_level_scaler - 5)
                        .attr("font-size", "10")
                        .attr("font-family", "Arial")
                        .text(function(){return Math.floor((max - min) / 1000) + "k"})
    
    for(i=0; i<maxloc-minloc; i+=200){
            if (i > 0 && i < maxloc && i % 2000 == 0){
                var tl = svg.append("text")
                            .attr("x", (minloc + i)/coefficient + 5)
                            .attr("y", y_level_scaler - 5)
                            .attr("font-size", "10")
                            .attr("font-family", "Arial")
                            .text(function(){return i / 1000 + "k"});
                var scaler = svg.append("line")
                            .attr("x1", (minloc + i)/coefficient + x_shifted)
                            .attr("y1", y_level_scaler - 3)
                            .attr("x2", (minloc + i)/coefficient + x_shifted)
                            .attr("y2", y_level_scaler + 3)
                            .attr("stroke", "black");
            } else {
                var scaler = svg.append("line")
                            .attr("x1", (minloc + i)/coefficient + x_shifted)
                            .attr("y1", y_level_scaler)
                            .attr("x2", (minloc + i)/coefficient + x_shifted)
                            .attr("y2", y_level_scaler + 2)
                            .attr("stroke", "black");
            };
    };
    
    // Draw backline for contigues
    var scaler = svg.append("line")
                    .attr("x1", minloc/coefficient + x_shifted)
                    .attr("y1", y_level_contig)
                    .attr("x2", maxloc/coefficient + x_shifted)
                    .attr("y2", y_level_contig)
                    .attr("stroke", "black")
                    .attr("id", "backline");
    
    // Draw contigues & address SNPs
    $.each(df, function(key, value){
        var points = "";
        if(value.strand == 1){
            points = (value.startloc / coefficient + x_shifted) + "," + (y_level_contig + 2.5) + " " + (value.stoploc / coefficient - 10  + x_shifted) + "," + (y_level_contig + 2.5) + " " + (value.stoploc / coefficient + x_shifted) + "," + (y_level_contig) + " " + (value.stoploc / coefficient - 10 + x_shifted) + "," + (y_level_contig - 2.5) + " " +  (value.startloc / coefficient + x_shifted) + "," + (y_level_contig - 2.5)
        }
        else{
            points = (value.startloc / coefficient + x_shifted) + "," + (y_level_contig) + " " + (value.startloc / coefficient + 10 + x_shifted) + "," + (y_level_contig - 2.5) + " " + (value.stoploc / coefficient + x_shifted) + "," + (y_level_contig - 2.5) + " " + (value.stoploc / coefficient + x_shifted) + "," + (y_level_contig + 2.5) + " " + (value.startloc / coefficient + 10 + x_shifted) + "," + (y_level_contig + 2.5)
        };
        var contig = svg.append("polygon")
                        .attr("points", points)
                        .style("fill", "darkgoldenrod")
                        .attr("id", key);
        
        $.each(value.snpsloc, function(varid, loc){
            var var_pos = svg.append("line")
                                .attr("x1", loc/coefficient + x_shifted)
                                .attr("y1", y_level_contig - 10)
                                .attr("x2", loc/coefficient + x_shifted)
                                .attr("y2", y_level_contig + 10)
                                .attr("stroke", "darkorchid")
                                .attr("class", "pos_mark")
                                .attr("id", varid);
            $("#" + varid + "_snp").click(function(){
                $(".pos_mark").not(this).hide();
                $("#" + varid).show();
                var Ref = $("#refflow");
                Ref.scrollLeft(loc/coefficient + x_shifted - 292);
            });
        })
    })
};

function drtree(df){
    var coefficient_x = 250;
    var coefficient_y = 32;
    var x_shifted = 30;
    var y_shifted = 60;
    var longest = 1.8;
    
    var svg = d3.select("#tree")
                .append("svg")
                .attr("height", 500)
                .attr("width", 500)
                .attr("id", "tree_svg");
    
    $.each(df.nodes, function(key, value){
        var x_branch = svg.append("line") 
                            .attr("x1", value.xcoord * coefficient_x + x_shifted)
                            .attr("y1", value.ycoord * coefficient_y + y_shifted)
                            .attr("x2", (value.xcoord - value.branch_length) * coefficient_x + x_shifted)
                            .attr("y2", value.ycoord * coefficient_y + y_shifted)
                            .attr("stroke", "black");
        if (value.is_Leaf != "1"){
            var y_branch = svg.append("line") 
                            .attr("x1", value.xcoord * coefficient_x + x_shifted)
                            .attr("y1", value.descendent_ycoord[0] * coefficient_y + y_shifted)
                            .attr("x2", value.xcoord * coefficient_x + x_shifted)
                            .attr("y2", value.descendent_ycoord[1] * coefficient_y + y_shifted)
                            .attr("stroke", "black");
        }
        else {
            var extension = svg.append("line") 
                                .attr("x1", value.xcoord * coefficient_x + x_shifted)
                                .attr("y1", value.ycoord * coefficient_y + y_shifted)
                                .attr("x2", longest * coefficient_x + x_shifted)
                                .attr("y2", value.ycoord * coefficient_y + y_shifted)
                                .attr("stroke", "red");
        };
        
    })
};

function dr_snp(df){
    // Draw gidlists
    var trs_for_gid = "";
    var trs = "<tr>";
    $.each(df, function(key, value){
        trs += "<td class='pos nuc_align " + value.snp_site + " " + value.strand + "' id='" + key + "_snp'>" + value.pos + "</td><td></td>";
    });
    trs += "</tr>";
    $.each(genome_id_list, function(i, d){
        trs += "<tr clas='" + d + "'>"
        trs_for_gid += "<tr class='" + d + "'><td>" + d + "</td></tr>"
        if (d == "gid_ref"){
            $.each(df, function(key, value){
                // apply reference SNPs
                trs += "<td class='gid_ref nuc'>" + value[d][0] + "</td>";
                trs += "<td class='gid_ref aa'>" + value[d][1] + "</td>";
            })            
        } else {
            $.each(df, function(key, value){
                trs += "<td class='gid_sample nuc " + value.snp_site + "'>" + value[d][0] + "</td>";
                trs += "<td class='gid_sample aa'>" + value[d][1] + "</td>";
            })
        };
        trs += "</tr>";
    });
    $("#gidnames").append(trs_for_gid);
    $("#DF").append(trs);
};