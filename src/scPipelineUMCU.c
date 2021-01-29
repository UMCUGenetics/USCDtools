#include <R.h>
#include <R_ext/Boolean.h>
#include <R_ext/Rdynload.h>
#include <Rdefines.h>
#include <htslib/sam.h>

/* -------------------------------------------------------------------------
 * Function: count_reads_for_range
 * ------------------------------------------------------------------------- */

SEXP
count_reads_for_range (SEXP input_file_sexp, SEXP region_sexp)
{
  bam_hdr_t  *bam_header        = NULL;
  htsFile    *bam_input_stream  = NULL;
  hts_idx_t  *bam_index         = NULL;
  hts_itr_t  *iterator          = NULL;
  const char *input_file        = NULL;
  const char *region            = NULL;
  SEXP       output             = NULL;
  int        *p_output          = NULL;
  int        reads_in_region    = 0;
  int        state              = 0;

  region     = CHAR(STRING_ELT(region_sexp, 0));
  input_file = CHAR(STRING_ELT(input_file_sexp, 0));

  bam_input_stream = hts_open (input_file, "r");
  if (! bam_input_stream)
    return Rf_ScalarLogical(FALSE);

  bam_index = sam_index_load (bam_input_stream, input_file);
  if (! bam_index)
    return Rf_ScalarLogical(FALSE);

  /* Read/write the SAM/BAM/CRAM header.
   * -------------------------------------------------------------------- */
  bam_header = sam_hdr_read (bam_input_stream);
  if (! bam_header)
    {
      hts_close (bam_input_stream);  bam_input_stream = NULL;
      hts_idx_destroy (bam_index);   bam_index = NULL;

      return Rf_ScalarLogical(FALSE);
    }

  /* Count the reads in the region.
   * -------------------------------------------------------------------- */
  iterator = sam_itr_querys (bam_index, bam_header, region);
  if (iterator != 0)
    {
      bam1_t *alignment = bam_init1 ();
      while ((state = sam_itr_next (bam_input_stream,
                                    iterator,
                                    alignment)) >= 0)
        {
	  reads_in_region++;
    	}

      bam_destroy1 (alignment); alignment = NULL;
    }

  sam_itr_destroy (iterator);    iterator = NULL;
  hts_idx_destroy (bam_index);   bam_index = NULL;
  bam_hdr_destroy (bam_header);  bam_header = NULL;
  hts_close (bam_input_stream);  bam_input_stream = NULL;

  PROTECT (output = NEW_INTEGER(1));
  p_output    = INTEGER_POINTER(output);
  p_output[0] = reads_in_region;
  UNPROTECT (1);

  return output;
}


/* -------------------------------------------------------------------------
 * Register the functions to R.
 * ------------------------------------------------------------------------- */

R_CallMethodDef callMethods[]  = {
  { "count_reads_for_range", (DL_FUNC)&count_reads_for_range, 2 },
  { NULL,                    NULL,                            0 }
};

void
R_init_read_counter (DllInfo *info)
{
  R_registerRoutines (info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols (info, FALSE);
}
