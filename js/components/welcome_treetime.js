import React from  'react'
import ReactDOM from 'react-dom'
import Header from './header.js'
var request = require('superagent');
import { Panel, PanelGroup, Button, Grid, Row, Col, FormControl, FormGroup, ControlLabel , Checkbox, Table, OverlayTrigger, Tooltip } from "react-bootstrap";
import GTR from './gtr.js'

var PanelText = React.createClass({
    render: function(){
        return (
            <div className="page_container">
                <div id="welcome_description">

                <h4>Features</h4>
                    <ul>
                        <li>Approximate maximum-likelihood time tree inference</li>
                        <li>Inference of GTR models</li>
                        <li>Rerooting to obtain best root-to-tip regression</li>
                        <li>Coalesent priors</li>
                        <li>Auto-correlated molecular clocks</li>
                    </ul>

                TreeTime source code is available on <a href='https://github.com/neherlab/treetime'>Github</a>.

                </div>

            </div>
        );
    }
});

var PanelFiles = React.createClass({

    onBuildTreeSelected: function(evt){
        var chk = this.props.TreeTimeConfig.build_tree != false && this.props.TreeTimeConfig.build_tree != null;
        var btree = !chk;
        this.props.setTreeTimeConfig({"build_tree":btree})
    },

    render: function(){
        const csv_tooltip = (
            <Tooltip className="csv_tooltip">
                Upload comma-separated file with sampling dates and additional metadata.
                The format of the file should be as follows:
                    <div className='spacer'/>
                    <table>
                    <tr>
                        <th>
                        name,date,location,...
                        </th>
                    </tr>
                    <tr>
                        <td>
                        A/Oregon/15/2009,2009.4,Oregon,... <br/>
                        A/New_York/182/2000,2000.1,New_York,...
                        </td>
                    </tr>
                    </table>
                    <div className='spacer'/>
                <ul className="scriptFont">
                <li>A single header row with column names is required. "name" and "date" are required fields.</li>

                <li>The first column should contain the <strong>name</strong> of the sequences.</li>

                <li>The first column with "date" in the name is interpreted as tip dates.
                Admissible formats are: (i) numeric date as "2015.7", (ii) date string, e.g. "YYYY-MM-DD".
                </li>

                <li>Other columns may contain arbitrary data of any format and have arbitrary names.
                The server will try to parse these columns and show them in the results section.</li>
                </ul>

            </Tooltip>
        );

        return (
            <div>
            <Panel collapsible defaultExpanded header="Upload data" className="panel-treetime" id="welcome_panel_files">
            <Grid id="welcome_upload_grid">

                {/*Upload Tree file*/}
                <Row className="grid-treetime-row">
                    <Col xs={6} md={4} id="welcome_col_upload_tree" className="grid-treetime-col-right" >
                        <span className="btn btn-primary btn-file btn-file-treetime" id="btn-1">
                            Newick
                            <input type="file" onChange={this.props.uploadTreeFile}/>
                        </span>
                        {this.props.appState.tree_filename}
                        <Checkbox
                            checked={this.props.TreeTimeConfig.build_tree}
                            onChange={this.onBuildTreeSelected}>
                            Build tree
                        </Checkbox>
                    </Col>
                </Row>

                {/*Upload Fasta file*/}
                <Row className="grid-treetime-row">
                    <Col  xs={6} md={4} className="grid-treetime-col-right">
                        <span className="btn btn-primary btn-file btn-treetime btn-file-treetime">
                            Fasta
                            <input type="file" onChange={this.props.uploadAlnFile} />
                        </span>
                        {this.props.appState.aln_filename}
                    </Col>
                </Row>

                {/*Upload Metadata file*/}
                <Row className="grid-treetime-row">
                    <Col xs={6} md={4} className="grid-treetime-col-right">
                        <OverlayTrigger placement="bottom" overlay={csv_tooltip}>
                        <span className="btn btn-primary btn-file btn-treetime btn-file-treetime">
                                CSV
                                <input type="file" onChange={this.props.uploadMetaFile} />
                        </span>
                        </OverlayTrigger>
                        {this.props.appState.meta_filename}
                    </Col>
                </Row>
            </Grid>
            </Panel>
            </div>
        );
    }
});


var PanelExamples = React.createClass({
    render: function(){
        return (
            <div>
            <Panel collapsible defaultCollapsed header="Example datasets" className="panel-treetime" id="welcome_panel_examples">
                            <Table striped condensed hove>
                            <thead>
                              <tr>
                                <th>Species</th>
                                <th>Genomic region</th>
                                <th>Length</th>
                                <th>Size</th>
                                <th>Date range</th>
                                <th></th>
                              </tr>
                            </thead>
                            <tbody>

                            {/*Influenza-H3N2-NA 20 seq example*/}
                            <tr>
                              <th>Influenza H3N2</th>
                              <th>NA</th>
                              <th>1409</th>
                              <th>20</th>
                              <th>2000-2013</th>
                              <th><Button bsStyle="primary" className="btn-treetime" onClick={this.props.onExample_H3N2_NA_20}>Load</Button></th>
                            </tr>

                            {/*Influenza-H3N2-HA 100 seq example*/}
                            <tr>
                              <th>Influenza H3N2</th>
                              <th>HA</th>
                              <th>1701</th>
                              <th>100</th>
                              <th>2011-2013</th>
                              <th><Button bsStyle="primary" className="btn-treetime" onClick={this.props.onExample_H3N2_HA_100}>Load</Button></th>
                            </tr>

                            {/*HIV RT 200 seq  example*/}
                            <tr>
                              <th>HIV subtype B</th>
                              <th>RT</th>
                              <th>1320</th>
                              <th>186</th>
                              <th>1978-2016</th>
                              <th><Button bsStyle="primary" className="btn-treetime" onClick={this.props.onExample_HIV_RT_200}>Load</Button></th>
                            </tr>

                            {/*HIV p17 200 seq  example*/}
                            <tr>
                              <th>HIV subtype B</th>
                              <th>p17</th>
                              <th>435</th>
                              <th>183</th>
                              <th>1978-2016</th>
                              <th><Button bsStyle="primary" className="btn-treetime" onClick={this.props.onExample_HIV_p17_200}>Load</Button></th>
                            </tr>

                            {/*Zika virus example*/}
                            <tr>
                              <th>Zika</th>
                              <th>Full genome</th>
                              <th>10617</th>
                              <th>65</th>
                              <th>2013-2016</th>
                              <th><Button bsStyle="primary" className="btn-treetime" onClick={this.props.onExample_zika_65}>Load</Button></th>
                            </tr>

                            {/*Ebola virus example*/}
                            <tr>
                              <th>Ebola</th>
                              <th>Full genome</th>
                              <th>19006</th>
                              <th>362</th>
                              <th>2014-2016</th>
                              <th><Button bsStyle="primary" className="btn-treetime" onClick={this.props.onExample_ebola}>Load</Button></th>
                            </tr>
                            </tbody>
                            </Table>

                        </Panel>

            </div>
        );
    }

});

var PanelConfig = React.createClass({

    getInitialState(){
        return ({
                available_gtrs:[]
            }
        );
    },

    componentWillMount: function() {
        console.log("Will mount: ");
        console.log(this.props.TreeTimeConfig.available_gtrs)
    },

    componentDidMount: function() {
        var dics = this.props.TreeTimeConfig;
        console.log("Did mount: ");
        console.log(dics)
    },

    componentWillUpdate(nextProps, nextState) {
        console.log("Will Update: ");
        console.log(nextProps)

        var new_gtrs = nextProps.TreeTimeConfig.available_gtrs;
        if (new_gtrs.length != this.state.available_gtrs.length){
            this.setState ({available_gtrs:nextProps.TreeTimeConfig.available_gtrs})
        }
    },

    onTreeTimeRoot : function (e){
        var chk = this.props.TreeTimeConfig.root != false && this.props.TreeTimeConfig.root != null;
        var new_root = chk ? false : 'best';
        this.props.setTreeTimeConfig({"root": new_root});
    },

    onTreeTimePoly : function(e){
        var chk = this.props.TreeTimeConfig.polytomies != false && this.props.TreeTimeConfig.polytomies != null;
        var new_poly = !chk;
        this.props.setTreeTimeConfig({"polytomies": new_poly});
    },

    onMuSelected : function(e){
        var chk = this.props.TreeTimeConfig.slope != false && this.props.TreeTimeConfig.slope != null;
        var use_slope = !chk;
        this.props.setTreeTimeConfig({"slope":use_slope})
    },

    onMuChanged : function (e){
        var val = e.target.value;
        this.props.setTreeTimeConfig({"slope_value":val})
    },

    onCoalescentPriorSelected : function(e){
        var chk = this.props.TreeTimeConfig.use_coalescent_prior != false &&
            this.props.TreeTimeConfig.use_coalescent_prior != null;
        var use_cp = !chk;
        this.props.setTreeTimeConfig({"use_coalescent_prior": use_cp})
    },
    onCoalescentPriorChanged : function(e){
        var val = e.target.value;
        this.props.setTreeTimeConfig({"coalescent_prior_value":val})
    },

    onRelaxClockSelected : function(e){
        var chk = this.props.TreeTimeConfig.use_relaxed_clock != false && this.props.TreeTimeConfig.use_relaxed_clock != null;
        var use_clock = !chk;
        this.props.setTreeTimeConfig({"use_relaxed_clock":use_clock})
    },

    onRelaxClockSlack : function(e){
        var val = e.target.value;
        var coupling = this.getRealxedClockCoupling()
        // NOTE have to update the whole relaxed clock dict
        this.props.setTreeTimeConfig({"relaxed_clock":{"slack": val, "coupling":coupling}})
    },

    onRelaxClockCoupling : function(e){
        var val = e.target.value;
        var slack = this.getRealxedClockSlack()
        // NOTE have to update the whole relaxed clock dict
        this.props.setTreeTimeConfig({"relaxed_clock":{"slack": slack, "coupling":val}})
    },

    getRealxedClockCoupling: function(){
        return this.props.TreeTimeConfig.relaxed_clock ?
                this.props.TreeTimeConfig.relaxed_clock.coupling :
                0;
    },

    getRealxedClockSlack: function(){
        return this.props.TreeTimeConfig.relaxed_clock ?
                this.props.TreeTimeConfig.relaxed_clock.slack :
                0;
    },

    elementShowStyle: function (show){
        return show ? {"display":"inline-block"} : {"display":"none"}
    },

    onGtrChange : function(e){

    },

    render: function(){
        const reroot_tooltip = (
            <Tooltip id="tooltip">
                Re-root tree to optimal root-to-tip regression.
            </Tooltip>
        );

        return (

            <div>
            <Panel collapsible defaultCollapsed header="Advanced configuration" className="panel-treetime" id="welcome_panel_config">

                {/*Reroot to best root*/}
                <Row>
                <Checkbox
                    checked={this.props.TreeTimeConfig.root != false && this.props.TreeTimeConfig.root != null}
                    onChange={this.onTreeTimeRoot}>
                    <OverlayTrigger placement="top" overlay={reroot_tooltip}>
                    <div>Optimize tree root position</div>
                    </OverlayTrigger>
                </Checkbox>
                </Row>

                {/*Polytomies resolution*/}
                <Row>
                <Checkbox
                    checked={this.props.TreeTimeConfig.polytomies != false && this.props.TreeTimeConfig.polytomies != null}
                    onChange={this.onTreeTimePoly}>
                    Resolve polytomies using temporal constraints
                </Checkbox>
                </Row>

                {/*GTR model*/}
                <Row>
                <FormGroup>
                    <GTR AppState={this.props.TreeTimeConfig} setTreeAncConfig={this.props.setTreeTimeConfig} setGtrState={this.props.setGtrState}/>
                </FormGroup>
                </Row>

                {/*Fix substitution rate*/}
                <Row>
                <FormGroup>
                    <Checkbox
                        onChange={this.onMuSelected}
                        checked={this.props.TreeTimeConfig.slope != false && this.props.TreeTimeConfig.slope != null}>
                        Fix substitution rate
                    </Checkbox>
                    <div style={this.elementShowStyle(this.props.TreeTimeConfig.slope)}>
                        <FormControl
                            type="number"
                            step="1e-4"
                            maxlength="5"
                            min="1e-9"
                            max="1e-2"
                            disabled={false}
                            style={{"display":"inline-block"}}
                            onChange={this.onMuChanged}
                            value={this.props.TreeTimeConfig.slope_value}/>
                        <span style={{"display":"inline-block"}}>(#/year)</span>
                    </div>
                </FormGroup>
                </Row>

                {/*Use coalescent prior*/}
                <Row>
                <FormGroup>
                    <Checkbox
                        checked={this.props.TreeTimeConfig.use_coalescent_prior != false && this.props.TreeTimeConfig.use_coalescent_prior != null}
                        onChange={this.onCoalescentPriorSelected}>
                        Use coalescent prior
                    </Checkbox>
                    <div style={this.elementShowStyle(this.props.TreeTimeConfig.use_coalescent_prior)}>
                        <FormControl
                            type="number"
                            step="1e-3"
                            maxlength="5"
                            min="0"
                            max="1e-1"
                            disabled={false}
                            style={{"display":"inline-block"}}
                            onChange={this.onCoalescentPriorChanged}
                            value={this.props.TreeTimeConfig.coalescent_prior_value}/>
                        <span style={{"display":"inline-block"}}>(Hamming distance)</span>
                    </div>
                </FormGroup>
                </Row>

                {/*Relaxed molecular clock*/}
                <Row>
                <FormGroup>
                    <Checkbox
                        onChange={this.onRelaxClockSelected}
                        checked={this.props.TreeTimeConfig.use_relaxed_clock != false && this.props.TreeTimeConfig.use_relaxed_clock != null}>
                        Relax molecular clock
                    </Checkbox>
                    <div style={this.elementShowStyle(this.props.TreeTimeConfig.use_relaxed_clock)}>
                        <div style={{"verticalAlign":"center", "lineHeight":"30px", "width":"100%"}}>
                            <span  style={{"display":"inline-block", "float":"left"}}>Slack: &alpha; =</span>
                             <FormControl
                                style={{"height":"30px", "display":"inline-block", "float":"left"}}
                                type="number"
                                min="0"
                                max="0.1"
                                step="1e-3"
                                onChange={this.onRelaxClockSlack}
                                value={this.getRealxedClockSlack()}>
                            </FormControl>
                            <span  style={{"display":"inline-block", "float":"left"}}>Coupling: &beta; =</span>
                            <FormControl
                                style={{"height":"30px", "display":"inline-block", "float":"left"}}
                                type="number"
                                min="0"
                                max="0.1"
                                step="1e-3"
                                value={this.getRealxedClockCoupling()}
                                onChange={this.onRelaxClockCoupling}>
                            </FormControl>
                        </div>
                    </div>
                </FormGroup>
                </Row>

                <Row>
                Additional features such as skyline inference are available in the treetime package but not accessible via the web interface.
                </Row>
            </Panel>
            </div>
        );
    }

});

var WelcomeTreeTimePage = React.createClass({

    getInitialState : function() {
      return {
        UID: null,
        // labels and status of the files uploads
        tree_file:false,
        tree_filename:"Select tree file",
        aln_file:false,
        aln_filename:"Select alignment file",
        meta_file:false,
        meta_filename:"Select meta data file",
        // treetime configuration
        TreeTimeConfig: {
            'available_gtrs':{},
            'gtr':null
        }
      };
    },

    componentDidMount: function(){
        var parentNode = this.getDOMNode().parentNode;
        // set user id from the HTML template
        this.state.UID = parentNode.attributes.userid.value;
        this.state.TreeTimeConfig = cfg;
        //this.setAppState({'UID':parentNode.attributes.userid.value});
        this.setState({'treetime_cfg':cfg})
    },

    setAppState : function(partialState){
        this.setState(partialState);
    },

    setTreeTimeConfig: function(cfg){
        var new_config = this.state.TreeTimeConfig
        for (var key in cfg) {
            new_config[key] = cfg[key];
            if (key == "build_tree" ){
                if (cfg[key]){
                    this.setState({
                    tree_filename: "Will be built from alignment",
                    tree_file: false,
                })
            }else{
                this.setState({
                    tree_filename: "Select tree file",
                    tree_file: false,
                })
            }
        }
        }
        this.setState({TreeTimeConfig:new_config})
        console.log("TreeTime config changed. New config: ");
        console.log(this.state.TreeTimeConfig)
    },

    // Files uploads section
    uploadTreeFile :function(evt){

        this.setTreeTimeConfig({"build_tree":false})

        if (evt.target.files.length == 0){
            console.log("Resetting treefile")
            this.setState({
                tree_filename:"Select tree file",
                tree_file:false,
            });
            return;
        }
        console.log("Uploading tree file...");
        var formData = new FormData();
        formData.append('treefile', evt.target.files[0]);
        //for (var key in evt.target.files) {
        //    formData.append(key, files[key]);
        //}


        this.setAppState({tree_file:true});
        request.post('/upload/treetime/' + this.state.UID + '/file')
            .send(formData)
            .end(this.onUploadTreeFile);
    },

    onUploadTreeFile: function(err, res){
        if (err){
            this.setAppState({tree_file:false, tree_filename: "Error uploading file"});
            alert("Tree file upload error. Please, try once more.")
            return;
        }

        this.setState({
            tree_filename:JSON.parse(res.text).TreeFile,
            tree_file:true
        })
    },


    uploadAlnFile :function(evt){

        if (evt.target.files.length == 0){
            // console.log("Resetting treefile")
            // this.setState({tree_filename:"No file chosen"})
            return;
        }
        //console.log("Uploading alignment file...");
        var formData = new FormData();
        formData.append('alnfile', evt.target.files[0]);
        //for (var key in evt.target.files) {
        //    formData.append(key, files[key]);
        //}
        this.setAppState({aln_file:true});
        request.post('/upload/treetime/' + this.state.UID + '/file')
            .send(formData)
            .end(this.onUploadAlnFile);
    },

    onUploadAlnFile: function(err, res){
        if (err){
            this.setAppState({aln_file:false});
            alert("Alignment file upload error. Please, try once more.")
            return;
        }
        this.setState({aln_filename:JSON.parse(res.text).AlnFile, aln_file:true})

    },

    uploadMetaFile :function(evt){

        if (evt.target.files.length == 0){
            // console.log("Resetting treefile")
            // this.setState({tree_filename:"No file chosen"})
            return;
        }
        //console.log("Uploading metadata file...");
        var formData = new FormData();
        formData.append('metafile', evt.target.files[0]);
        //for (var key in evt.target.files) {
        //    formData.append(key, files[key]);
        //}

        this.setAppState({meta_file:true});
        request.post('/upload/treetime/' + this.state.UID + '/file')
            .send(formData)
            .end(this.onUploadMetaFile);
    },

    onUploadMetaFile: function(err, res){
        if (err){
            this.setAppState({meta_file:false});
            alert("Meta data file upload error. Please, try once more.")
            return;
        }

        this.setState({meta_filename:JSON.parse(res.text).MetaFile, meta_file:true})
    },

    //Examples section
    onExample_H3N2_NA_20 : function(){
        this.runExample("H3N2_NA_20");
    },

    onExample_H3N2_HA_100 : function(){
        console.log("Running example for H3N2_100")
        this.runExample("H3N2_HA_100")
    },

    onExample_HIV_RT_200: function(){
        this.runExample("HIV_RT_200")
    },

    onExample_HIV_p17_200: function(){
        this.runExample("HIV_p17_200")
    },

    onExample_H3N2_NA_500 : function(){
        this.runExample("H3N2_NA_500");
    },

    onExample_zika_65 : function(){
        this.runExample("zika_65");
    },

    onExample_ebola : function(){
        this.runExample("ebola");
    },

    runExample : function(example){
        console.log("run example requested:  " + example)
        this.setTreeTimeConfig({"build_tree":false})
        request.post("/treetime/" + this.state.UID + "/example")
          .set('Content-Type', 'application/json')
          .send({example:example})
          .end(this.onExampleUpload);
    },

    onExampleUpload : function(err, res){

        console.log("Example files uploaded")
        this.onUploadAlnFile(err, res)
        this.onUploadMetaFile(err, res)
        this.onUploadTreeFile(err, res)

    },


    onRunTreeTime: function(){
        //console.log("APP:: RUN button pressed");
        if ((!this.state.tree_file & !this.state.TreeTimeConfig.build_tree) || ! this.state.aln_file || !this.state.meta_file){
          var msg = "Cannot proceed with TreeTime: one or more file not loaded.\n\n"
          if ((!this.state.tree_file & !this.state.TreeTimeConfig.build_tree)){
            msg += "Phylogenetic tree file is missing.\n\n";
          }
          if (!this.state.aln_file){
            msg += "Sequence alignment file is missing.\n\n";
          }
          if (!this.state.meta_file){
            msg += "Meta data file is missing.\n\n";
          }
          alert(msg);
          return;
        }
        request.post("/treetime/" + this.state.UID + "/run")
          .set('Content-Type', 'application/json')
          .send({config: this.state.TreeTimeConfig})
          .end(this.onRunTreeTimeResponse);
    },

    onRunTreeTimeResponse : function(err, res){

        //console.log("RUN RESPONSE");
        //console.log(res)
        if (err){
            alert("Cannot start TreeTime calculations. Server error.")
            console.log(err)
            return;
        }
        window.location.replace("/treetime/" + this.state.UID + "/progress");
    },

    setGtrState: function(key, param_name, param_value){
        console.log("Welcome page: setting GTR state")
        var cfg = this.state.TreeTimeConfig
        var gtr = cfg.available_gtrs[key]
        if (!gtr.params){
            alert("Cannot set GTR parameter: This GTR has no parameters.")
            return;
        }
        for (var i = 0; i < gtr.params.length; ++i){
            var param = gtr.params[i];
            if (param.name == param_name){
                this.state.TreeTimeConfig.available_gtrs[key].params[i].value = param_value
                this.forceUpdate()
                break;
            }
        }
    },

    render:function(){
        return (
            <div>
                <Header/>
                <div className="page_container">
                <div className='bigspacer'/>

                {/* Upload files*/}
                <PanelFiles
                    TreeTimeConfig={this.state.TreeTimeConfig}
                    setTreeTimeConfig={this.setTreeTimeConfig}
                    appState={this.state}
                    uploadTreeFile={this.uploadTreeFile}
                    uploadAlnFile={this.uploadAlnFile}
                    uploadMetaFile={this.uploadMetaFile}
                />

                {/* Choose predefined example dataset*/}
                <PanelExamples
                    onExample_H3N2_NA_20={this.onExample_H3N2_NA_20}
                    onExample_H3N2_HA_100={this.onExample_H3N2_HA_100}
                    onExample_HIV_RT_200={this.onExample_HIV_RT_200}
                    onExample_HIV_p17_200={this.onExample_HIV_p17_200}
                    onExample_zika_65={this.onExample_zika_65}
                    onExample_ebola={this.onExample_ebola}
                />

                {/* Advanced configuration*/}
                <PanelConfig
                    setGtrState={this.setGtrState}
                    TreeTimeConfig={this.state.TreeTimeConfig}
                    setTreeTimeConfig={this.setTreeTimeConfig}/>

                {/* Run treetime on the server*/}
                <Button bsStyle="primary" onClick={this.onRunTreeTime}>Run treetime</Button>
                </div>
            </div>
        );
    }
});

ReactDOM.render((
    <WelcomeTreeTimePage/>),
    document.getElementById('react'));

export default WelcomeTreeTimePage;
