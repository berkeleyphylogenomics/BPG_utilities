$(function() {

    //This should probably be available via the API, but this is quicker right now.
    // (If you are reading this, you can interpret that as an apology - sorry)
    // each field can have different operators and values, or an input box, I guess
    var filter_templates = [
        {
            "field": "accession",
            "title": "UniProt Accession",
            "operators": [ "exact" ],
            "type": "input"
        },
        {
            "field": "uniprot_identifier",
            "title": "UniProt Identifier",
            "operators": [ "exact" ],
            "type": "input"
        },
        {
            "field": "description",
            "title": "Description",
            "operators": [ "contains", "not_contains" ],
            "type": "input"
        },
        {
            "field": "taxon__id",
            "title": "Species",
            "operators": [ "exact" ],
            "type": "select",
            "options" : [
                { "name": "Human", "value": "9606" },
                { "name": "Escherichia coli", "value": "562" },
                { "name": "Arabidopsis pumila", "value": "74718" },
                { "name": "Saccharomyces cerevisiae", "value": "4932" },
                { "name": "Bacillus subtilis", "value": "1423" }
            ]
        },
        {
            "field": "in_swissprot_f",
            "title": "In SwissProt",
            "operators": [ "exact" ],
            "type": "select",
            "options" : [
                { "name": "Yes" , "value": "1" },
                { "name": "No" , "value": "0" }
            ]
        },
        {
            "field": "go_exp_evidence",
            "title": "GO Experimental Evidence",
            "operators": [ "exact" ],
            "type": "select",
            "options" : [
                { "name": "Yes" , "value": "1" },
                { "name": "No" , "value": "0" }
            ]
        },
        {
            "field": "go_annotations__go_term__accession",
            "title": "GO Accession (GO:0000011)",
            "operators": [ "exact" ],
            "type": "input"
        },
        {
            "field": "ec_annotations__ec__description",
            "title": "Enzyme Commission Description",
            "operators": [ "contains", "not_contains" ],
            "type": "input"
        },
        {
            "field": "ec_annotations__ec__class_number",
            "title": "Enzyme Commission Class Number",
            "operators": [ "exact" ],
            "type": "input"
        },
        {
            "field": "ec_annotations__ec__subclass_number",
            "title": "Enzyme Commission Sub Class Number",
            "operators": [ "exact" ],
            "type": "input"
        },
        {
            "field": "ec_annotations__ec__subsubclass_number",
            "title": "Enzyme Commission Sub Sub Class Number",
            "operators": [ "exact" ],
            "type": "input"
        },
        {
            "field": "ec_annotations__ec__enzyme_number",
            "title": "Enzyme Commission Enzyme Number",
            "operators": [ "exact" ],
            "type": "input"
        },
        {
            "field": "uniprotpfam__pfam__accession",
            "title": "Pfam Accession (PF08793)",
            "operators": [ "exact" ],
            "type": "input"
        },
        {
            "field": "uniprotpfam__pfam__name",
            "title": "Pfam Name (NB-ARC)",
            "operators": [ "exact" ],
            "type": "input"
        }
    ];


    // Takes a field hash, backbone cid and dumps out an HTML string object
    // why isn't this a template?  Because it was "simple"? maybe the conditional stuff
    function RenderField(field, cid) {
        var html = "<div class='new-field'>";

        html += "<input type='radio' name='field' value='" + field.field + "'>";
        html += "<input type='hidden' id='" + cid + "_" + field.field + "_title' name='title' value='" + field.title + "'>";

        html += "<label for='" + cid + "_" + field.field +"'>";
        html += field.title;
        html += "</label>";

        html += "<select id='" + cid + "_" + field.field + "_operator' name='operator'>";
        html += _.map(field.operators, function(operator) { return "<option>" + operator + "</option>"; }).join('');
        html += "</select>";

        switch (field.type) {
            case "select":
                html += "<select id='" + cid + "_" + field.field + "_value' name='value'>";
                html += _.map(field.options, function(option) { return "<option value='" + option.value + "'>" + option.name + "</option>"; }).join('');
                html += "</select>";
                break;
            case "input":
                html += "<input id='" + cid + "_" + field.field + "_value' name='value' >";
                break;
        }

        html += "</div>";
        return html;
    }


    // So, why backbone?  basically events and consistency.
    // doing standard encapsulated objects works really well from a data perspective,
    // but these objects really do the same thing as the view layer for the most part
    // being able to painlessly link up the events was what made me actually do it.
    // it isn't hard, but why reinvent the wheel?


    FilterClause = Backbone.View.extend({
        template: _.template($('#filter_clause_template').html()),
        tagname: 'div',
        initialize: function(arg) {
            // The name of the field in the REST API
            this.field    = arg.field,
            // The DB Operator (we just use the django ones, like exact, containes, lte, etc
            this.operator = arg.operator,
            // Comparison Value in the DB
            this.value    = arg.value;
            // pretty title text for the "humans" who don't like seeing db column names
            this.title    = arg.title;
            // We have to traverse up to delete ourselves
            this.parent = arg.parent;
        },
        // returns an html string suitible for DOM insertion
        render: function () {
            $(this.el).html(this.template(this));
            $(this.el).find("#remove_" + this.cid).click( $.proxy(this.remove, this) )
            return this;
        },

        // returns a json fragment
        serialize: function () {
            return {
                "field":    this.field,
                "operator": this.operator,
                "value":    this.value
            };
        },
        remove: function() {
            this.parent.remove_child(this);
        }
    });
    
    FilterSet = Backbone.View.extend({
        template: _.template($('#filter_set_template').html()),
        tagname: 'div',
        add_filter_set: function() {
            var operator = $(this.el).find(".new-filter-set [name=operator]:checked").val();
            this.children.push( new FilterSet({children: [], operator: operator, parent: this }) );
            return this.render();
        },
        add_filter_clause: function() {
            var field = $(this.el).find(".new-filter-clause input[name=field]:checked").val();
            var operator = $(this.el).find("#"+this.cid+"_"+field+"_operator").val();
            var value = $(this.el).find("#"+this.cid+"_"+field+"_value").val();
            var title = $(this.el).find("#"+this.cid+"_"+field+"_title").val();
            this.children.push( new FilterClause({field: field, operator: operator, value: value, title: title, parent: this }) );
            return this.render();
        },
        initialize: function(arg) {
            // These can be either other FilterSets, or they can be FilterClauses
            this.children = arg.children || [],
        
            // How the set is reduced.  Either AND or OR, at least to start
            this.operator = arg.operator;

            // We have to traverse up to delete ourselves
            this.parent = arg.parent;

        },
        render: function () {
            var self = this;
            $(this.el).html(
                this.template({
                    "filters": _.map( filter_templates, function(filter) { return RenderField(filter, self.cid) }).join(''),
                    "cid": this.cid
                })
            );
            _.each( this.children, function(child) { 
                $(self.el).find("#children_"+self.cid).append( child.render().el ); 
            });

            // Activate the event handlers (remember the earlier comment about how I'm using backbone because their event framework
            // It sometimes makes things much simpler, I promise.. really..
            $(this.el).find("#add_filter_set_" + this.cid).click( $.proxy(this.add_filter_set, this) )
            $(this.el).find("#add_filter_clause_" + this.cid).click( $.proxy(this.add_filter_clause, this) )
            $(this.el).find("#remove_" + this.cid).click( $.proxy(this.remove, this) )

            return this;
        },
        remove: function() {
            this.parent.remove_child(this);
        },
        remove_child: function(child) {
            idx = this.children.indexOf(child);
            this.children.splice(idx, 1);
            delete child;
            this.render();
        },

        // serializes all children for our server request
        serialize: function () {
            return {
                "operator": this.operator,
                "children": _.map(this.children, function(child) { return child.serialize(); } )
            };
        }
    });

    Query = Backbone.View.extend({
        el: "#query",

        initialize: function() {
            // A query is made up of chained filter sets.
            // I think.  This may be subject to change once I realize in what precise way I am an idiot.
            this.filter_sets = [];
            this.operator = "IN";

            // Is this just a special case of a filter set?  It really seems like it
            this.filter_sets.push( new FilterSet({  children : [ ], operator: "AND" }) );

            this.render();
        },
        render: function() {
            $(this.el).find("#query_description").empty();
            // There are only a max of two elements, and there isn't an easy reverse iterator, sorry
            this.filter_sets.reverse();
            _.each( this.filter_sets, function(filter, idx) {
                if (idx !== 0 ) {
                    $("#query_description").append("<p class='collect'>With an ortholog that has these properties</p>");
                }
                $( filter.render().el ).appendTo($("#query_description") );
            });
            this.filter_sets.reverse();
        },
        serialize: function() {
            var filters = _.map(this.filter_sets, function(filter_set) { return filter_set.serialize(); } );
            var json_filter = $.toJSON( filters );
            return json_filter;
        },
        events: {
            "click #add_orthologs": "add_orthologs"
        },
        add_orthologs: function() {
            this.filter_sets.push( new FilterSet({  children : [ ], operator: "AND" }) );
            this.filter_sets.reverse();
            this.render();
            $("#add_orthologs").hide();
            create_table({with_phogs: true});
        }
    });

    function make_sequence_link(acc_or_ident) {
        return "<a href='/phylofacts/sequence/UniProt/" + acc_or_ident + "' >" + acc_or_ident + "</a>";
    }

    function make_phog_link(acc) {
        return "<a href='/phog/" + acc + "' >" + acc + "</a>";
    }

    function create_table(params) {
        params = params || {};
        //TODO - replace this with something better
        if (params.with_phogs) {
            columns = [
                { "sTitle": "Accession", "mDataProp": "accession","fnRender": function( o, val ) { return make_sequence_link(o.aData.accession); } },
                { "sTitle": "Identifier", "mDataProp": "uniprot_identifier", "fnRender": function( o ) { return make_sequence_link(o.aData.uniprot_identifier); }  },
                { "sTitle": "Description", "mDataProp": "de" },
                { "sTitle": "PHOG", "mDataProp": "phog_accession", "fnRender": function( o, val ) { return make_phog_link(o.aData.phog_accession); } },
                { "sTitle": "PHOG's UniProt Accession", "mDataProp": "phog_uniprot_accession", "fnRender": function( o, val ) { return make_sequence_link(o.aData.phog_uniprot_accession); } },
                { "sTitle": "PHOG's UniProt Description", "mDataProp": "phog_uniprot_de" },
            ];
        } else {
            columns = [
                { "sTitle": "Accession", "mDataProp": "accession","fnRender": function( o, val ) { return make_sequence_link(o.aData.accession); } },
                { "sTitle": "Identifier", "mDataProp": "uniprot_identifier", "fnRender": function( o ) { return make_sequence_link(o.aData.uniprot_identifier); }  },
                { "sTitle": "Description", "mDataProp": "de" },
            ];
        }

        $("#results").empty().html("<table></table>");
        ResultsTable = $("#results table").dataTable({
            "bJQueryUI": true, 
            "bFilter": false,
            "bProcessing": true,
            "bServerSide": true,
            "sAjaxSource": '/api/sequence/',
            "aoColumns": columns,
            "fnServerParams": function ( aoData ) {
                aoData.push( { "name": "filters", "value": QueryInstance.serialize() } );
            }
        });
        table_created = 1;
    }

    QueryInstance = new Query;
    var table_created = 0;
    $("#refresh").button().click(function() { 
        if (!table_created) {
            create_table();
            table_created = 1;
        }
        ResultsTable.fnDraw() 
    });
    $("#add_orthologs").button();



});
