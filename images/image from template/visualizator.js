class myContainer
{
    constructor(playgroud) {
        this.type = "MyContainer";
        this.playgroud = playgroud;
    }
    setData(allData){
        this.allData = allData;
        console.log("All data is set");
        this.startDrawing();
    }
    startDrawing(){
        var self = this;
        // using self.allData = {   template: {},
        //                          draftContainer: [{}, {}, {links: [], params: {}}]}
        var traces = [];
        _.each(self.allData["draftContainer"], function (t) { 
            let p = t["params"];
            var trace = {
                x: t["links"],
                type: "histogram",
                name: p["label"],
                opacity: p["alpha"],
                marker: {
                    color: p["facecolor"],
                },
                histnorm: "probability density",
                histfunc: "count",
                xbins:{
                    size: self.allData["template"]["binSize"],
                },
            }
            traces.push(trace);
        });
        var layout = {
            barmode: "overlay",
            title: this.allData["template"]["title"], 
            xaxis: {title: this.allData["template"]["xlabel"]}, 
            yaxis: {title: this.allData["template"]["ylabel"]}
        };
        
        this.P = Plotly.newPlot(this.playgroud, traces, layout);
    }

}


