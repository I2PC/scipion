/*
 * File:        Transpose.js
 * Version:     1.0.0
 * CVS:         $Id$
 * Description: Transpose (rotate) a table about the diagonal axis, turning columns into rows and vice versa
 * Author:      Fargin Bastage
 * Created:     Sun Aut 7 11:23:29 CDT 2011
 * Modified:    $Date$ by $Author$
 * Language:    Javascript
 * License:     GPL v2 or BSD 3 point style
 * Project:     DataTables
 * Dependency:  JQuery, DataTables
 * Contact:     fbastage@yahoo.com
 * 
 * Copyright 2011 Fargin Bastage, all rights reserved.
 *
 * This source file is free software, under either the GPL v2 license or a
 * BSD style license, available at:
 *   http://datatables.net/license_gpl2
 *   http://datatables.net/license_bsd
 *
 */


(function($) {

/* fnTransposeVersion
 * 
 * arguments:
 *   oSettings, DataTables automatically puts table settings object as first param
 *   sVersion, optional version string to check current Transpose version against 
 *     (passing a non-string, like an int or float results in TypeError)
 *
 * returns:
 *   if sVersion is not provided, or is empty string, 0, or false, return sTransposeVersion
 *   else return true if sTransposeVersion is greater than or equal to sVersion (treated as series of numbers)
 *   else return false
 *
 * note: for comparison, version string is considered a series of integers
 *   only the first 3 numbers will be considered.
 *   1.0.0.1500 will be equal to 1.0.0.1
 *
 * only positive integers should be used. letters at the end of numbers will be gracefully ignored
 *   i.e. 1.-3.5 will have undefined behavior
 *   1.a.5  will be treated as 1.0.5 - as thus, it should be avoided to avert confusion
 *   1.3a.005complete  will be a valid version string, equivalent to 1.3.5
 *
 */
$.fn.dataTableExt.oApi.fnTransposeVersion = function ( oSettings, sVersion )
{
  if ( typeof oSettings == "undefined" || typeof sVersion == "undefined" || !sVersion) return Transpose.VERSION;
  
  // else 
  if (sVersion == Transpose.VERSION) return true;
  
  aVersion = sVersion.split('.');
  aTransposeVersion = Transpose.VERSION.split('.');
  
  // convert version numbers to int
  for (i = 0; i < aVersion.length; i++) aVersion[i] = parseInt(aVersion[i]);
  for (i = 0; i < aTransposeVersion.length; i++) aTransposeVersion[i] = parseInt(aTransposeVersion[i]);
  
  if (aTransposeVersion.length < 1) return true;
  if (aVersion.length < 1) return true;
  if (aTransposeVersion[0] > aVersion[0]) return true;
  else if (aTransposeVersion[0] < aVersion[0]) return false;
  
  if (aTransposeVersion.length < 2) return true;
  if (aVersion.length < 2) return true;
  if (aTransposeVersion[1] > aVersion[1]) return true;
  else if (aTransposeVersion[1] < aVersion[1]) return false;
  
  if (aTransposeVersion.length < 3) return true;
  if (aVersion.length < 3) return true;
  if (aTransposeVersion[2] > aVersion[2]) return true;
  else if (aTransposeVersion[2] < aVersion[2]) return false;
  
  // if we have gotten this far, version numbers are considered equal
  return true;
  
}

/* Transpose or untranspose the table 
 *
 * (if bTranspose is not provided, use true)
 */
$.fn.dataTableExt.oApi.fnTranspose = function ( oSettings, bTranspose )
{
  if (typeof bTranspose == "undefined") bTranspose = true;
  oSettings.oInstance.oPluginTranspose.fnTranspose(bTranspose);
}

/* Get / set transpose state (without triggering redraw)
 *
 */
$.fn.dataTableExt.oApi.fnTransposeState= function ( oSettings, bTranspose )
{
  if (typeof bTranspose != "undefined")
    oSettings.oInstance.oPluginTranspose.s.bTranspose = bTranspose;
  
  return oSettings.oInstance.oPluginTranspose.s.bTranspose;
}

/* 
 * Add functionality feature to DataTables
 *
 * aoInitComplete: turn display of DataTable to "none" 
 * aoDrawCallback: perform redraw of transpose table on every table draw
 *
 */
$.fn.dataTableExt.aoFeatures.push( {
  "fnInit": function( oSettings ) {
    var oTable = oSettings.oInstance;
    if ( typeof oTable.oPluginTranspose == 'undefined' ) {
      var opts = typeof oSettings.oInit.oPluginTranspose != 'undefined' ? 
        oSettings.oInit.oPluginTranspose : {};
      oTable.oPluginTranspose = new Transpose( oSettings.oInstance, opts );
    } else {
      oTable.oApi._fnLog( oSettings, 1, "Transpose attempted to initialise twice. Ignoring second." );
    }
        
    oSettings.aoDrawCallback.push( {
      "sName": "Transpose",
      "fn": function () {
        oTable.oPluginTranspose.fnTranspose();
      }      
    } );  
    
    return null; /* No node to insert */
  },
  "cFeature": "Z",
  "sFeature": "Transpose"
} );


/** 
 * Transpose provides rotating a DataTable about it's diagonal axis so columsn become rows and vice versa
 * @class Transpose
 * @constructor
 * @param {object} DataTables object
 * @param {object} Transpose options
 */
Transpose = function( oTable, oOpts )
{
  /* Santiy check that we are a new instance */
  if ( !this.CLASS || this.CLASS != "Transpose" )
  {
    alert( "Warning: Transpose must be initialised with the keyword 'new'" );
  }
  
  if ( typeof oOpts == 'undefined' )
  {
    oOpts = {};
  }
  
  
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   * Public class variables
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  
  /**
   * @namespace Settings object which contains customisable information for Transpose instance
   */
  this.s = {
    /**
     * DataTables settings object
     *  @property dt
     *  @type     Object
     *  @default  null
     */
    "dt": null,
    
    /**
     * Initialisation object used for this instance
     *  @property init
     *  @type     object
     *  @default  {}
     */
    "init": oOpts,
    
    /**
     * Whether or not to transpose the table
     *  @property bTranspose
     *  @type     int
     *  @default  true
     */
    "bTranspose": true,
    
    /**
     * Whether or not to transpose the table
     *  @property bTranspose
     *  @type     int
     *  @default  true
     */
    "sTransposeId": null
  };
   
 
  
  /* Constructor logic */
  this.s.dt = oTable.fnSettings();
  this._fnConstruct();
  
  return this;
};




Transpose.prototype = {
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   * Public methods
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  
  "fnTranspose": function (bTranspose)
  {
    if (typeof bTranspose != "undefined") this.s.bTranspose = bTranspose && true;
    
    if (this.s.bTranspose) this._fnShowTranspose();
    else this._fnHideTranspose();
    
  },
  
  
  /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
   * Private methods (they are of course public in JS, but recommended as private)
   * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
  
  /**
   * Constructor logic
   *  @method  _fnConstruct
   *  @returns void
   *  @private 
   */
  "_fnConstruct": function ()
  {
    var that = this;
    
    /* Transpose id */
    if ( typeof this.s.init.sTransposeId != 'undefined' )
    {
      this.s.sTransposeId = this.s.init.sTransposeId;
    }
    else if ( typeof this.s.dt.oInit.sTransposeId != 'undefined' )
    {
      this.s.sTransposeId = this.s.dt.oInit.sTransposeId;
    }    
    else if ( typeof this.s.dt.sTableId != 'undefined' )
    {
      this.s.sTransposeId = this.s.dt.sTableId + "_transpose";
    } 
    
    
    /* Transpose state */
    if ( typeof this.s.init.bTranspose != 'undefined' )
    {
      this.s.bTranspose = this.s.init.bTranspose;
    }
    else if ( typeof this.s.dt.oInit.bTranspose != 'undefined' )
    {
      this.s.bTranspose = this.s.dt.oInit.bTranspose;
    }
    else 
    {
      this.s.bTranspose = true;
    }
  },
  
  /**
   * Transpose DataTable about diagonal axis. (create rotated copy of table and hide original table)
   *  Calling this function when transposed table is already shown will re-draw the table
   *  @method  _fnShowTranspose
   *  @returns void
   *  @private 
   */
  "_fnShowTranspose": function ()
  {
    
    var oSettings = this.s.dt;
    var sTableId = this.s.dt.sTableId;
    var nTable = this.s.dt.nTable;
    $(nTable).hide();  // hide the data table
    $('.dataTables_scrollHeadInner').hide();  // hide scroll head (when sScrollX is used)
    $('.dataTables_scrollFootInner').hide();  // hide scroll foot (when sScrollX is used)

    // get any previous transposed tables, to destroy after creating a new one
    var sTransposeId = this.s.sTransposeId;
    var old_transpose = $('#'+sTransposeId);

    // create new transpose table from original table
    var transpose = $(nTable).clone(true);
    var $transpose = $(transpose);

    $transpose.addClass('DTTZ_Transposed');
    $transpose.attr('id', sTransposeId);
    $transpose.find('tr').remove();  // remove all rows
    $tbody = $($transpose.find('tbody')[0]);

    // rebuild table row by row by reading a column header and all column data from the original table
    for (i = 1; i <= oSettings.aoColumns.length; i++) {  // nodes are 1-indexed in :nth-child
        // start a new row
        var $nTr = $('<tr>');
        $nTr.addClass(oSettings.asStripeClasses[i % oSettings.asStripeClasses.length]);

        // when sScrollX is used, the header is placed into a separate div. get the header and place into column 0
        if ($('.dataTables_scrollHeadInner').length) {
          $('.dataTables_scrollHeadInner thead th:nth-child('+i+')').each( function() {
          	console.log("header2")
        	console.log($(this))
        	console.log("finheader2")
            $nTr.append($(this).clone(true)); 
          });
        } else {
          $('thead tr th:nth-child('+i+')', nTable).each( function() { 
            	console.log("header3")
            	console.log($(this))
            	console.log("finheader3")  
                    	console.log("vdengana44444444444n")
            $nTr.append($(this).clone(true)); 
          	console.log($nTr) 
          });
        }

        // get columns from original table, add to new row
        $('tbody tr td:nth-child('+i+')', nTable).each( function() { 
            $nTr.append($(this).clone(true)); 
        });

        $tbody.append($nTr[0]);
    }
    


    if (old_transpose.length) old_transpose.remove();  // remove any old versions of transpose table, from old draws
    $transpose.append($tbody[0]);

    $(nTable).after($transpose[0]);  // add the newly created table
    $transpose.show();  //show, because it inherited from hidden table.    
    
  },
  
  
  /**
   * Untranspose DataTable. (destroy transpoed copy of table and show original table)
   *  @method  _fnHideTranspose
   *  @returns void
   *  @private 
   */
  "_fnHideTranspose": function ()
  {
     $('#'+this.s.sTransposeId).remove();
     $('#'+this.s.dt.sTableId).show();
     $('.dataTables_scrollHeadInner').show();  // hide scroll head (when sScrollX is used)
     $('.dataTables_scrollFootInner').show();  // hide scroll foot (when sScrollX is used)

  }
  
  
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Constants
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**
 * Name of this class
 *  @constant CLASS
 *  @type     String
 *  @default  Transpose
 */
Transpose.prototype.CLASS = "Transpose";


/**
 * ColReorder version
 *  @constant  VERSION
 *  @type      String
 *  @default   As code
 */
Transpose.VERSION = "1.0.0";
Transpose.prototype.VERSION = Transpose.VERSION;


})(jQuery); 
