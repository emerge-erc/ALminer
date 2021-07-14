Query keywords
==============

Below is a table of possible keywords that can be used to query the ALMA
Science Archive using the `alminer.keysearch`_ function:

.. list-table::
   :widths: 25 25 50
   :header-rows: 1
   :align: left

   * - ALMA query keyword
     - Type
     - Description
   * - access_format
     - char(9)
     - Content format of the data         
   * - access_url
     - char(72*)
     - URL to download the data
   * - antenna_arrays
     - char(660*)
     - Blank-separated list of Pad:Antenna pairs, i.e., A109:DV09 J504:DV02 J505:DV05 for antennas DV09, DV02 and DV05 sitting on pads A109, J504, and J505, respectively.
   * - asdm_uid
     - char(32*)
     - UID of the ASDM containing this Field.                                                             
   * - authors
     - char(4000*)
     - Full list of first author and all co-authors                                                       
   * - band_list
     - char(30*)
     - Space delimited list of bands
   * - bib_reference
     - char(30*)
     - Bibliography code
   * - data_rights
     - char(11)
     - Access to data.
   * - dataproduct_type
     - char(5*)
     - type of product
   * - facility_name
     - char(3)
     - telescope name
   * - first_author
     - char(256*)
     - The first author as provided by `telbib.eso.org`_.
   * - frequency_support
     - char(4000*)
     - All frequency ranges used by the field                                                             
   * - group_ous_uid
     - char(64*)
     - Group OUS ID
   * - instrument_name
     - char(4)
     - instrument name
   * - is_mosaic
     - char(1)
     - Flag to indicate if this ASDM represents a mosaic or not.                                          
   * - lastModified   
     - char(*)
     - Time stamp of last modification of the metadata                                                    
   * - member_ous_uid 
     - char(64*)   |
     - Member OUS ID
   * - o_ucd    
     - char(35)    |
     - UCD describing the observable axis (pixel values)                                                  
   * - obs_collection 
     - char(4)
     - short name for the data collection
   * - obs_creator_name    | 
     - char(256*)  |
     - case-insensitive partial match over the full PI name. Wildcards can be used                        
   * - obs_id   
     - char(64*)   |
     - internal dataset identifier
   * - obs_publisher_did
     - char(33*)
     - publisher dataset identifier
   * - obs_release_date
     - char(*)
     - timestamp of date the data becomes publicly available                                              
   * - obs_title
     - char(256*)
     - Case-insensitive search over the project title                                                     
   * - pol_states
     - char(64*)
     - polarization states present in the data                                                            
   * - proposal_abstract
     - char(4000*)
     - Text search on the proposal abstract. Only abstracts will be returned which contain the given text. The search is case-insensitive.                                     
   * - proposal_authors
     - char(2000*)
     - Full name of CoIs.
   * - proposal_id
     - char(64*)
     - Identifier of proposal to which NO observation belongs.                                            
   * - pub_abstract
     - char(4000*)
     - Case insensitive text search through the abstract of the publication.                              
   * - pub_title
     - char(256*)
     - Case insensitive search through the title of the publication.                                      
   * - qa2_passed
     - char(1)
     - Quality Assessment 2 status: does the Member / Group OUS fulfil the PI's requirements?
   * - s_region 
     - char(*)
     - region bounded by observation
   * - scan_intent
     - char(256*)
     - Scan intent list for the observed field.                                                           
   * - schedblock_name
     - char(128*)
     - Name of the Scheduling Block used as a template for executing the ASDM containing this Field.      
   * - science_keyword
     - char(200*)
     - Chosen by the PI in the observing tool at the time of proposal submission. For an overview, see `Appendix D of the ALMA Proposer's Guide`_. For a precise list, see a `this table of science keywords`_.
   * - science_observation
     - char(1)
     - Flag to indicate whether this is a science observation.                                            
   * - scientific_category
     - char(200*)
     - Chosen by the PI in the observing tool at the time of proposal submission. For an overview, see `Appendix D of the ALMA Proposer's Guide`_. For a precise list, see `this table of scientific categories`_.
   * - target_name
     - char(256*)
     - name of intended target
   * - type
     - char(16*)
     - Type flags.

.. _alminer.keysearch: ../pages/api.rst#alminer.keysearch
.. _telbib.eso.org: http://telbib.eso.org
.. _Appendix D of the ALMA Proposer's Guide: https://almascience.eso.org/proposing/proposers-guide#section-63
.. _this table of scientific categories: ../pages/scientific_categories.rst
.. _this table of science keywords: ../pages/science_keywords.rst