/*
 
 FUNCTION TO READ GFF/GTF ANNOTATION FILES AND FILTER POSITIONS AND VARIANTS DEFINED BY THE USER:
 
 FIRST PART: READ THE FILE GFF/GTF, INCLUDE ALL ROWS IN A STRUCTURE AND FILTER ACCORDING THE OVERLAPPED CRITERIA TRANSCRIPTS. 
 SECOND PART: ONCE THE GFF ROWS ARE INCLUDED AND FILTERED ACCORDING A OVERLAPPED TRANSCRIPTS CRITERIA, DETECT POSITIONS AND VARIANTS

 Sebastian E. Ramos-Onsins
*/


#include "usegff.h"

static char tripletsN[64][3] =
{
	{"111"},	{"112"},	{"114"},	{"113"},	{"121"},	{"122"},	{"124"},	{"123"},
	{"141"},	{"142"},	{"144"},	{"143"},	{"131"},	{"132"},	{"134"},	{"133"},
	{"211"},	{"212"},	{"214"},	{"213"},	{"221"},	{"222"},	{"224"},	{"223"},
	{"241"},	{"242"},	{"244"},	{"243"},	{"231"},	{"232"},	{"234"},	{"233"},
	{"411"},	{"412"},	{"414"},	{"413"},	{"421"},	{"422"},	{"424"},	{"423"},
	{"441"},	{"442"},	{"444"},	{"443"},	{"431"},	{"432"},	{"434"},	{"433"},
	{"311"},	{"312"},	{"314"},	{"313"},	{"321"},	{"322"},	{"324"},	{"323"},
	{"341"},	{"342"},	{"344"},	{"343"},	{"331"},	{"332"},	{"334"},	{"333"},
};
	
int use_gff(char *name_fileinputgff,char *subset_positions,char *genetic_code,
			double *matrix_sizepos,int n_samp,long int n_site,char *DNA_matr,
			double *matrix_segrpos,FILE *file_output,int mainargc,int include_unknown,
			char *criteria_transcripts, int type_output, long int *nmhits, long int *mhitbp,
			int outgroup_presence, int nsamoutg)
{
	
	FILE *file_gff,*file_gff2;
	char *row,*f,cstrand[1],cframe[2],aaseq[1],aaput[1];
	long int i=0;
	int nn,q,qq,stop,countpath,xx,countpath3;
	long int jj,l,j,n,m,nrows,jo,hf1,hf2;
	char fields[9][1024];	
	struct valuesgff *fieldsgff,*fieldsgff2;
	char *seqid/*, *fileid*/;
	double *cmat,*cmatnc,*cmatsil;
	char *cframe_pos;
	long int ii,ii2,k;
	long int startframe	= 0;
	long int endframe	= 0;
	long int endtrp		= 0;
	char *cod3n			= 0;
	char *cod3put		= 0;
	long int start,end;
	double cod3f,cod3ft[3];
	int neffsam;
	
	int **read_frame;
	
	int gg;
	char *pst;
	int *vector_erase_overlapped;
	
	long int ibeg,iend,nrows2;
	int overlap_nf;
	double rest=0.;
	double rest2=0.;
	long int ntransc;
	char **matrix_coding;
	char transcript[256],let[1];
	char name_fileinputgff2[1024];
	char *ff,*ff1;
	
	int countframe;
	long int *longtr;
	long int lotr;
    long int ncountrow;
	
    int n_fold[64][3];/*count the number of degenerate times for each position of the codon*/
    
    cstrand[0]=0;
	cframe[0]=0;
	cframe[1]=0;
	aaseq[0]=0;
	aaput[0]=0;

    /*calculate the degeneration of each triplet and position:*/
    if(function_do_nfold_triplets(n_fold,genetic_code,tripletsN) == 0)
    {
        fprintf(file_output,"\nError: It is not possible to define degenerated positions. use_gff.c\n");
        return(0);
    }

    /*
    FIRST PART: READ THE FILE GFF/GTF, INCLUDE ALL ROWS IN A STRUCTURE AND FILTER ACCORDING THE OVERLAPPED CRITERIA TRANSCRIPTS.
	 
	 */
	/*read gff file: read name_fileinputfoldergff, if it does not work, read name_fileinputfolderGFF*/
	/*name_fileinputgff/name_fileinputGFF is only for displaying the name on the screen.*/
	/*fields to read in GFF-format: filename, noread(source), feature, start, end, noread(score), strand, frame, [gene_id. transcript_id, Parent, noread(rest)]*/

	if((file_gff = fopen (name_fileinputgff,"r")) == 0)
	{
		printf("\n  It is not possible to open the file %s",name_fileinputgff);
		return 0; /*error*/
	}

	if(file_gff) {
		/*init*/
		if( !(f = (char *) malloc( BUFSIZ ))) { /**< TODO: Check que BUFSIZ es el de stdio */
			fprintf(file_output,"\nError: memory not reallocated. use_gff.1 \n");
			return 0; /*error*/
		}
		if(!(row = (char *)malloc(1024*sizeof(char)))) {
			fprintf(file_output,"\nError: memory not reallocated. use_gff.2 \n");
			return 0; /*error*/
		}
		if(!(fieldsgff = (struct valuesgff *)calloc(1,sizeof(struct valuesgff)))) {
			fprintf(file_output,"\nError: memory not reallocated. use_gff.3 \n");
			return 0; /*error*/
		}
		setbuf(file_gff, f);
		/*read rows*/
		nrows = 0;
        ncountrow = 0;
		while(1)
		{
            ncountrow += 1;
			row[i=0] = '\0';
			fgets(row, 1024*sizeof(char), file_gff);
			if(row[i] == '\0' && feof(file_gff)) break;
			/*i=0;*/
			/*while((row[i] = fgetc(file_gff)) != 0 && row[i] != 10 && row[i] != 13 && feof(file_gff)!= 1) {*/
			/*	i++;*/
			/*	if(i >= 1024) while((j = fgetc(file_gff)) != 0  && j != 10 && j != 13 && (feof(file_gff)) != 1); */
			/*}*/
			/*row[i] = 0;*/
			
			i=0;
			while(/*row[i] == 32 || */row[i] == '\t') {
				if(row[i] == 10 || row[i] == 13 || row[i] == 0) break;
				i++;
				if(i >= 1024) break;
			}
			if(row[i] == '#') continue;
			if(i >= 1024) continue;
			/*include fields in variables*/
			j = k = 0;
			while(row[i] != 10 && row[i] != 13 && row[i] != 0) {
				while(/*row[i] == 32 || */row[i] == '\t') {
					if(row[i] == 10 || row[i] == 13 || row[i] == 0) break;
					i++;
					if(i >= 1024) break;
				}
				k=0;
				while(/*row[i] != 32 && */row[i] != '\t' && row[i] != 10 && row[i] != 13 && row[i] != 0) {
					fields[j][k] = row[i];
					k++;
					if(k >= 256) break;
					i++;
					if(i >= 1024) break;
				}
				fields[j][k] = '\0';
				j++;
				if(j==9) {
					break;
				}
			}
			if(i < 16) {
				continue;
			}
			
			/*add fields in struct*/
			for(n=0;n<j;n++) {
				switch(n) {
					case 0:
						fieldsgff[nrows].filename[0] = '\0';
						strncat(fieldsgff[nrows].filename,fields[n],256*sizeof(char));
						break;
					case 1:
						fieldsgff[nrows].source[0] = '\0';
						strncat(fieldsgff[nrows].source,fields[n],256*sizeof(char));
						break;
					case 2:
						fieldsgff[nrows].feature[0] = '\0';
						strncat(fieldsgff[nrows].feature,fields[n],256*sizeof(char));
						break;
					case 3:
						fieldsgff[nrows].start = atol(fields[n]);
						break;
					case 4:
						fieldsgff[nrows].end = atol(fields[n]);
						break;
					case 5:
						fieldsgff[nrows].score[0] = '\0';
						strncat(fieldsgff[nrows].score,fields[n],256*sizeof(char));
						break;
					case 6:
						fieldsgff[nrows].strand[0] = fields[n][0];
						break;
					case 7:
						if(fields[n][0] == '1') {
							fieldsgff[nrows].frame[0] = '2';
						}
						else {
							if(fields[n][0] == '2') {fieldsgff[nrows].frame[0] = '1';
							}
							else {
								fieldsgff[nrows].frame[0] = fields[n][0];
							}
						}
						break;
					case 8:
                        fieldsgff[nrows].gene_id[0] = '\0';
                        pst = strstr(fields[n],"ID=");
                        if(pst) {
                            while(*pst != '=') pst++;
                            pst++;
                            gg = 0;
                            while(*pst != ';' && *pst != '\0') {
                                fieldsgff[nrows].gene_id[gg] = *pst;
                                pst++;
                                gg++;
                                if(gg==255) break;
                                
                            }
                            fieldsgff[nrows].gene_id[gg] = '\0';
                        }
						fieldsgff[nrows].transcript_id[0] = '\0';
						pst = strstr(fields[n],"Parent=");
						if(pst) {
							while(*pst != '=') pst++;
							pst++;
							gg = 0;
							while(*pst != ';' && *pst != '\0') {
								fieldsgff[nrows].transcript_id[gg] = *pst;
								pst++;
								gg++;
								if(gg==255) break;
								
							}
							fieldsgff[nrows].transcript_id[gg] = '\0';
						}
						fieldsgff[nrows].seqname[0] = '\0';
						pst = strstr(fields[n],"Target=");
						if(pst) {
							while(*pst != '=') pst++;
							pst++;
							gg = 0;
							while(*pst != ';' && *pst != '\0') {
								fieldsgff[nrows].seqname[gg] = *pst;
								pst++;
								gg++;
								if(gg==255) break;
								
							}
							fieldsgff[nrows].seqname[gg] = '\0';
						}
                        if(fieldsgff[nrows].gene_id[0] == '\0') {
                            pst = strstr(fields[n],"gene_id=");
                            if(pst) {
                                while(*pst != '=') pst++;
                                pst++;
                                gg = 0;
                                while(*pst != ';' && *pst != '\0') {
                                    fieldsgff[nrows].gene_id[gg] = *pst;
                                    pst++;
                                    gg++;
                                    if(gg==255) break;
                                    
                                }
                                fieldsgff[nrows].gene_id[gg] = '\0';
                            }
                        }
                        if(fieldsgff[nrows].gene_id[0] == '\0') {
                            pst = strstr(fields[n],"gene_id ");
                            if(pst) {
                                while(*pst != 32) pst++;
                                pst++;
                                gg = 0;
                                while(*pst != ';' && *pst != '\0') {
                                    fieldsgff[nrows].gene_id[gg] = *pst;
                                    pst++;
                                    gg++;
                                    if(gg==255) break;
                                    
                                }
                                fieldsgff[nrows].gene_id[gg] = '\0';
                            }
                        }
                        if(fieldsgff[nrows].gene_id[0] == '\0') {
                            pst = strstr(fields[n],"name ");
                            if(pst) {
                                while(*pst != 32) pst++;
                                pst++;
                                gg = 0;
                                while(*pst != ';' && *pst != '\0') {
                                    fieldsgff[nrows].gene_id[gg] = *pst;
                                    pst++;
                                    gg++;
                                    if(gg==255) break;
                                    
                                }
                                fieldsgff[nrows].gene_id[gg] = '\0';
                            }
                        }
                        if(fieldsgff[nrows].transcript_id[0] == '\0') {
                            pst = strstr(fields[n],"transcript_id ");
                            if(pst) {
                                while(*pst != 32) pst++;
                                pst++;
                                gg = 0;
                                while(*pst != ';' && *pst != '\0') {
                                    fieldsgff[nrows].transcript_id[gg] = *pst;
                                    pst++;
                                    gg++;
                                    if(gg==255) break;
                                    
                                }
                                fieldsgff[nrows].transcript_id[gg] = '\0';
                            }
                        }
                        if(fieldsgff[nrows].transcript_id[0] == '\0') {
                            pst = strstr(fields[n],"transcript_id=");
                            if(pst) {
                                while(*pst != '=') pst++;
                                pst++;
                                gg = 0;
                                while(*pst != ';' && *pst != '\0') {
                                    fieldsgff[nrows].transcript_id[gg] = *pst;
                                    pst++;
                                    gg++;
                                    if(gg==255) break;
                                    
                                }
                                fieldsgff[nrows].transcript_id[gg] = '\0';
                            }
                        }
                        if(fieldsgff[nrows].transcript_id[0] == '\0') {
                            pst = strstr(fields[n],"proteinId ");
                            if(pst) {
                                while(*pst != 32) pst++;
                                pst++;
                                gg = 0;
                                while(*pst != ';' && *pst != '\0') {
                                    fieldsgff[nrows].transcript_id[gg] = *pst;
                                    pst++;
                                    gg++;
                                    if(gg==255) break;
                                    
                                }
                                fieldsgff[nrows].transcript_id[gg] = '\0';
                            }
                        }
						
						if(strcmp(fieldsgff[nrows].seqname,fieldsgff[nrows].gene_id) != 0) {
							if(fieldsgff[nrows].gene_id[0] != '\0') strcpy(fieldsgff[nrows].seqname,fieldsgff[nrows].gene_id);
							else strcpy(fieldsgff[nrows].gene_id,fieldsgff[nrows].seqname);
						}
						
						break;
				}
			}
			/*filtering ...*/
			if(fieldsgff[nrows].start < (long int)1 &&
               fieldsgff[nrows].strand != '\0')  {
				if(type_output == 0 || type_output == 10) {
					if(file_output) fprintf(file_output,"GFF file error: start (%ld) is lower than 1. Row %ld not analyzed. ",fieldsgff[nrows].start,ncountrow);
				}
				continue;
			}
			if(fieldsgff[nrows].end > (long int)n_site &&
               fieldsgff[nrows].strand != '\0')  {
				if(type_output == 0 || type_output == 10) {
					if(file_output) fprintf(file_output,"GFF file error: end (%ld) is larger than number of total sites (%ld). Row %ld not analyzed.\n ",fieldsgff[nrows].end,n_site,ncountrow);
				}
				continue;
			}
			if(fieldsgff[nrows].start > fieldsgff[nrows].end && fieldsgff[nrows].strand != '\0') {
				if(type_output == 0 || type_output == 10) {
					if(file_output) fprintf(file_output,"GFF file error: start (%ld) is larger than end (%ld). Row %ld not analyzed. ",fieldsgff[nrows].start,fieldsgff[nrows].end,ncountrow);
				}
				continue;
			}
			/*internally starts from 0 to n_site-1*/
			fieldsgff[nrows].start -= (long int)1;
			fieldsgff[nrows].end -= (long int)1;
			/*we assume by default starting from strand '+' and frame 0*/
			/*minimum five fields: seqname,whatever(dot),feature,start,end.*/
			if(j < 9) continue;
			/**/
			if(strcmp(subset_positions,"synonymous") == 0 || strcmp(subset_positions,"nonsynonymous") == 0 || strcmp(subset_positions,"silent") == 0 || strcmp(subset_positions,"0-fold") == 0 || strcmp(subset_positions,"2-fold") == 0 || strcmp(subset_positions,"3-fold") == 0 || strcmp(subset_positions,"4-fold") == 0) {
				if(strcmp(fieldsgff[nrows].feature,"CDS") == 0) {
					if(fieldsgff[nrows].strand[0] != '\0' && fieldsgff[nrows].strand[0] != '+' && 
					   fieldsgff[nrows].strand[0] != '-'  && fieldsgff[nrows].strand[0] != '.') {
						fieldsgff[nrows].strand[0] = 32;
					}
					if(fieldsgff[nrows].frame[0] != '\0' && fieldsgff[nrows].frame[0] != '1' && 
					   fieldsgff[nrows].frame[0] != '2'  && fieldsgff[nrows].frame[0] != '0' && 
					   fieldsgff[nrows].frame[0] != '.') {
						fieldsgff[nrows].frame[0] = 32;
					}
				}
			}
			/**/
			/*if row accepted*/
			nrows += 1;
			if(!(fieldsgff = (struct valuesgff *)realloc(fieldsgff,(nrows+1)*sizeof(struct valuesgff)))) {
				fprintf(file_output,"\nError: memory not reallocated. use_gff.3 \n");
				free(row);
				free(f);
				free(fieldsgff);
				fclose(file_gff);
				return 0; /*error*/
			}		
			if(feof(file_gff)) break;
		}
		free(row);
		free(f);
		fclose(file_gff);
		
		for(n=0;n<9;n++) {
			switch(n) {
				case 0:
					fieldsgff[nrows].filename[0] = '\0';
					strncat(fieldsgff[nrows].filename,name_fileinputgff,256*sizeof(char));
					break;
				case 1:
					fieldsgff[nrows].source[0]='\0';
					break;
				case 2:
					fieldsgff[nrows].feature[0] = '\0';
					strncat(fieldsgff[nrows].feature,"all",256*sizeof(char));
					break;
				case 3:
					fieldsgff[nrows].start = (long int)0;
					break;
				case 4:
					fieldsgff[nrows].end = (long int)(n_site-1);
					break;
				case 5:
					fieldsgff[nrows].score[0]='\0';
					break;
				case 6:
					fieldsgff[nrows].strand[0] = 32;
					break;
				case 7:
					fieldsgff[nrows].frame[0] = 32;
					break;
				case 8:
					fieldsgff[nrows].seqname[0] = '\0';
					fieldsgff[nrows].gene_id[0] = '\0';
					fieldsgff[nrows].transcript_id[0] = '\0';
					break;
			}
		}
		nrows += 1;
		
		if(!(read_frame = (int **)calloc(1,sizeof(int *)))) {
			fprintf(file_output,"\nError: memory not reallocated. use_gff.34 \n");
			free(fieldsgff);
			return 0; /*error*/
		}
		if(!(read_frame[0] = (int *)calloc(256,sizeof(int)))) {
			fprintf(file_output,"\nError: memory not reallocated. use_gff.341 \n");
			free(read_frame);
			free(fieldsgff);
			return 0; /*error*/
		}
		
		/*set to zero all values in matrix_sizepos*/
		if(strcmp(subset_positions,"noncoding") == 0) for(ii=0;ii<n_site;ii++) matrix_sizepos[ii] = (double)1;
		else for(ii=0;ii<n_site;ii++) matrix_sizepos[ii] = (double)0;
		/*init a current matrix_sizepos for each seqname/gene_id*/
		if((cmat = (double *)calloc((unsigned long)n_site,sizeof(double))) == 0) {
			fprintf(file_output,"\nError: memory not reallocated. use_gff.4 \n");
			return 0; /*error*/
		}
		if((cmatnc = (double *)calloc((unsigned long)n_site,sizeof(double))) == 0) {
			fprintf(file_output,"\nError: memory not reallocated. use_gff.4b \n");
			return 0; /*error*/
		}
		if((cmatsil = (double *)calloc((unsigned long)n_site,sizeof(double))) == 0) {
			fprintf(file_output,"\nError: memory not reallocated. use_gff.4c \n");
			return 0; /*error*/
		}
		for(ii=0;ii<n_site;ii++) cmatsil[ii] = (double)1;
		/*init a 3*n_samp matrix for syn/nsyn triplets*/
		if(strcmp(subset_positions,"synonymous") == 0 || strcmp(subset_positions,"nonsynonymous") == 0 || strcmp(subset_positions,"silent") == 0 || strcmp(subset_positions,"0-fold") == 0 || strcmp(subset_positions,"2-fold") == 0 || strcmp(subset_positions,"3-fold") == 0 || strcmp(subset_positions,"4-fold") == 0) {
			if((cod3n = (char *)calloc((unsigned long)(3*(n_samp+1)),sizeof(char))) == 0) {
				fprintf(file_output,"\nError: memory not reallocated. use_gff.5 \n");
				return 0; /*error*/
			}
			if((cod3put = (char *)calloc((unsigned long)3,sizeof(char))) == 0) {
				fprintf(file_output,"\nError: memory not reallocated. use_gff.5 \n");
				return 0; /*error*/
			}
		}
		
		/**********************************************  BEGIN TO CHECK  *******************************************************/
		/**********************************************  BEGIN TO CHECK  *******************************************************/
		/**********************************************  BEGIN TO CHECK  *******************************************************/
		/**********************************************  BEGIN TO CHECK  *******************************************************/
		/**********************************************  BEGIN TO CHECK  *******************************************************/
		/**********************************************  BEGIN TO CHECK  *******************************************************/
		/**********************************************  BEGIN TO CHECK  *******************************************************/
		/**********************************************  BEGIN TO CHECK  *******************************************************/
		
		/*
		 once all rows are counted and kept, look at each gene_id transcripts. Take the "criteria_transcripts" (max or min).
		 - given the criteria (max,min) keep only the transcripts following the criteria. if overlapped, create a new one merging all.
		 - in case overlapping regions out of coding frame, eliminate ALSO this region in matrix_sizepos (0 in those positions, it is not noncoding).
		 */
		
		/*create a new struct and copy all rows*/
		if(!(fieldsgff2 = (struct valuesgff *)calloc(nrows,sizeof(struct valuesgff)))) {
			fprintf(file_output,"\nError: memory not reallocated. use_gff.3b \n");
			return 0; /*error*/
		}
		/*sort the struct by gene_id(first), start and transcript_id*/
		qsort(fieldsgff,nrows/*-1*/,sizeof(struct valuesgff),comp_gene_id);
		
		/* keep exons from the same gene_id/seqname according overlapping criteria*/
		/*first define all '.' values*/
		for(n=0;n<nrows-1;) {
			seqid = fieldsgff[n].seqname;
			if(seqid[0]=='\0') seqid = fieldsgff[n].gene_id;
			if(fieldsgff[n].gene_id[0] =='\0') memcpy(fieldsgff[n].gene_id,fieldsgff[n].seqname,sizeof(char)*256);
			m = n+1;
			while(m<nrows && strcmp(fieldsgff[m].gene_id,seqid) == 0) {m++;} /*from n to m: all rows with the same gene_id*/
			if(n==nrows) break;
			
			/*Positive strand: first define all CDS rows that have '.' by its real value in frame and strand.*/
			qsort(fieldsgff+n,m-n,sizeof(struct valuesgff),comp_start_id); /*sort within gene_id by start*/
			qsort(fieldsgff+n,m-n,sizeof(struct valuesgff),comp_trcpt_id); /*sort within gene_id by transcript_id*/
			
			/*considering a given gene_id*/
			for(j=n;j<m;j++) {
				while(j<m && strcmp(fieldsgff[j].feature,"CDS") != 0) {j++;} /*look for CDS rows*/
				if(j==m) break;
				seqid = fieldsgff[j].transcript_id;
				cstrand[0] = fieldsgff[j].strand[0];
				cframe[0]  = fieldsgff[j].frame[0];
				if(cstrand[0] == '.' && j==n) { /*j row is undefined*/
					/*printf("\n Reading GTF file: CDS exons with undefined sense strand (gene_id = %s, transcript_id = %s). TRANSCRIPT NOT CONSIDERED.",fieldsgff[j].gene_id,fieldsgff[j].transcript_id);*/
					strcpy(fieldsgff[j].feature,"CDS_EXCLUDED\0");
					for(ii=j+1;ii<m;ii++) if(strcmp(fieldsgff[m].gene_id,seqid) == 0 && strcmp(fieldsgff[ii].feature,"CDS") == 0) {strcpy(fieldsgff[ii].feature,"CDS_EXCLUDED\0");}
					continue;
				}
				if(cframe[0] == '.' && cstrand[0] == '+' && j==n) {
					/*printf("\n Reading GTF file: CDS exons with undefined coding frames (gene_id = %s, transcript_id = %s). TRANSCRIPT NOT CONSIDERED.",fieldsgff[j].gene_id,fieldsgff[j].transcript_id);*/
					strcpy(fieldsgff[j].feature,"CDS_EXCLUDED\0");
					for(ii=j+1;ii<m;ii++) if(strcmp(fieldsgff[m].gene_id,seqid) == 0 && strcmp(fieldsgff[ii].feature,"CDS") == 0) {strcpy(fieldsgff[ii].feature,"CDS_EXCLUDED\0");}
					continue;
				}
				ii = j;
				
				k = j+1;
				while(k<m && strcmp(fieldsgff[k].transcript_id,seqid) == 0) {k++;} /*from j to k: rows with the same transcript_id*/				
				
				/*considering a given transcript_id*/
				for(l=j+1;l<k;l++) {
					while(l<k && strcmp(fieldsgff[l].feature,"CDS") != 0) {l++;} /*look for CDS rows*/
					if(l==k) break;
					if(cstrand[0] != fieldsgff[l].strand[0] && fieldsgff[l].strand[0] != '.') {
						/*printf("\n Reading GTF file: CDS exons with different sense strand (gene_id = %s, transcript_id = %s). TRANSCRIPT NOT CONSIDERED.",fieldsgff[j].gene_id,fieldsgff[j].transcript_id);*/
						strcpy(fieldsgff[j].feature,"CDS_EXCLUDED\0");
						for(ii=j+1;ii<k;ii++) {
							if(strcmp(fieldsgff[ii].transcript_id,seqid) == 0 && strcmp(fieldsgff[ii].feature,"CDS") == 0) {
								strcpy(fieldsgff[ii].feature,"CDS_EXCLUDED\0");
							}
						}
						break;
					}
					
					if(cstrand[0] != fieldsgff[l].strand[0] && fieldsgff[l].strand[0] == '.') {
						fieldsgff[l].strand[0] = cstrand[0];
					}
					if(fieldsgff[l].frame[0] == '.' && fieldsgff[l].strand[0] == '+') {
						rest = fmod((double)( fieldsgff[ii].end-fieldsgff[ii].start),3.0);
						rest = fmod((double)((fieldsgff[ii].frame[0]-48) + rest),3.0);       
						fieldsgff[l].frame[0] = 48 + (char)fmod((double)rest+1.0,3.0);
					}
					ii = l;
				}
			}
			
			/*negative strand, the frame is the opposite: starts in the position end of the exon located more on the left*/
			qsort(fieldsgff+n,m-n,sizeof(struct valuesgff),comp_end_id); /*sort within gene_id by end*/
			qsort(fieldsgff+n,m-n,sizeof(struct valuesgff),comp_trcpt_id); /*sort within gene_id by transcript_id*/
			
			for(j=m-1;j>=n;j--) {
				while(j>=n && strcmp(fieldsgff[j].feature,"CDS") != 0) {j--;} /*look for CDS rows*/
				if(j<n) break;
				seqid = fieldsgff[j].transcript_id;
				cstrand[0] = fieldsgff[j].strand[0];
				cframe[0]  = fieldsgff[j].frame[0];
				if(cstrand[0] == '.' && j==m-1) { /*j row is undefined*/
					/*printf("\n Reading GTF file: CDS exons with undefined sense strand (gene_id = %s, transcript_id = %s). TRANSCRIPT NOT CONSIDERED.",fieldsgff[j].gene_id,fieldsgff[j].transcript_id);*/
					strcpy(fieldsgff[j].feature,"CDS_EXCLUDED\0");
					for(ii=j-1;ii>=n;ii--)
						if(strcmp(fieldsgff[m].gene_id,seqid) == 0 && strcmp(fieldsgff[ii].feature,"CDS") == 0) {
							strcpy(fieldsgff[ii].feature,"CDS_EXCLUDED\0");
						}
					continue;
				}
				if(cframe[0] == '.' && cstrand[0] == '-' && j==m-1) {
					/*printf("\n Reading GTF file: CDS exons with undefined coding frames (gene_id = %s, transcript_id = %s). TRANSCRIPT NOT CONSIDERED.",fieldsgff[j].gene_id,fieldsgff[j].transcript_id);*/
					strcpy(fieldsgff[j].feature,"CDS_EXCLUDED\0");
					for(ii=j-1;ii>=n;ii--) 
						if(strcmp(fieldsgff[m].gene_id,seqid) == 0 && strcmp(fieldsgff[ii].feature,"CDS") == 0) {
							strcpy(fieldsgff[ii].feature,"CDS_EXCLUDED\0");
						}
					continue;
				}
				ii = j;
				
				k = j;
				while(k>n && strcmp(fieldsgff[k].transcript_id,seqid) == 0) {k--;} /*from j to k: rows with the same transcript_id*/				
				
				/*considering a given transcript_id*/
				for(l=k;l>=j;l--) {
					while(l>=j && strcmp(fieldsgff[l].feature,"CDS") != 0) {l--;} /*look for CDS rows*/
					if(l<j) break;
					if(cstrand[0] != fieldsgff[l].strand[0] && fieldsgff[l].strand[0] != '.') {
						/*printf("\n Reading GTF file: CDS exons with different sense strand (gene_id = %s, transcript_id = %s). TRANSCRIPT NOT CONSIDERED.",fieldsgff[j].gene_id,fieldsgff[j].transcript_id);*/
						strcpy(fieldsgff[j].feature,"CDS_EXCLUDED\0");
						for(ii=k;ii>=j;ii--) {
							if(strcmp(fieldsgff[ii].transcript_id,seqid) == 0 && strcmp(fieldsgff[ii].feature,"CDS") == 0) {
								strcpy(fieldsgff[ii].feature,"CDS_EXCLUDED\0");
							}
						}
						break;
					}
					
					if(cstrand[0] != fieldsgff[l].strand[0] && fieldsgff[l].strand[0] == '.') {
						fieldsgff[l].strand[0] = cstrand[0];
					}
					if(fieldsgff[l].frame[0] == '.' && fieldsgff[l].strand[0] == '-') {
						rest = fmod((double)( fieldsgff[ii].end-fieldsgff[ii].start),3.0);
						rest = fmod((double)((fieldsgff[ii].frame[0]-48) + rest),3.0); 
						fieldsgff[l].frame[0] = 48 + (char)fmod((double)rest+1.0,3.0);
					}
					ii = l;
				}
			}
			n = m;
		}	
		
		nrows2 = 0;
		for(n=0;n<nrows;) {
			while(n<nrows && strcmp(fieldsgff[n].feature,"CDS") != 0) {
				memcpy(fieldsgff2+nrows2,fieldsgff+n,sizeof(struct valuesgff)*1);/*include row in fieldsgff2*/
				nrows2++;
				/*
				if(nrows2 >= nrows) {
					if(!(fieldsgff2 = (struct valuesgff *)realloc(fieldsgff2,(nrows2+1)*sizeof(struct valuesgff)))) {
						fprintf(file_output,"\nError: memory not reallocated. use_gff.3b \n");
						return 0;
					}
				}
				*/
				n++;
			}
			if(n==nrows) break;
			
			seqid = fieldsgff[n].gene_id;
			m = n+1;
			while(m<nrows && strcmp(fieldsgff[m].gene_id,seqid) == 0) {m++;} /*from n to m: all rows with the same gene_id*/			
			/*if(m==nrows) break;*/
			/*considering a given gene_id*/
			
			/*how many transcripts has this gene_id?*/
			qsort(fieldsgff+n,m-n,sizeof(struct valuesgff),comp_start_id); /*sort within gene_id by start*/
			qsort(fieldsgff+n,m-n,sizeof(struct valuesgff),comp_trcpt_id); /*sort within gene_id by transcript_id*/
			j = n;
			while(strcmp(fieldsgff[j].feature,"CDS")!=0) {j++;}
			seqid = 0; /*fieldsgff[j].transcript_id;*/
			ntransc=0; 
			j = n;
			while(j<m) {
				while(j<m && strcmp(fieldsgff[j].feature,"CDS")!=0) {j++;}
				if(j<m && strcmp(fieldsgff[j].feature,"CDS") == 0 && ((seqid == 0) || strcmp(fieldsgff[j].transcript_id,seqid) != 0)) {
					ntransc++;
					seqid = fieldsgff[j].transcript_id;
					j++;
				}
				while(j<m && strcmp(fieldsgff[j].transcript_id,seqid) == 0 && strcmp(fieldsgff[j].feature,"CDS") == 0) {j++;}
			}
			
			/*define char matrix of coding regions x transcripts of the gene_id located in n*/
			start = n_site;
			end = 0;
			for(j=n;j<m;j++) {
				if(strcmp(fieldsgff[j].feature,"CDS")==0) {
					if(start > fieldsgff[j].start) start = fieldsgff[j].start;
					if(end   < fieldsgff[j].end)   end   = fieldsgff[j].end;
				}
			}
			if(!(matrix_coding = (char **)calloc(ntransc+1,sizeof(char *)))) {
				fprintf(file_output,"\nError: memory not reallocated. use_gff.23b \n");
				return 0; /*error*/
			}
			j=n;
			for(i=0;i<ntransc;i++) {
				if(!(matrix_coding[i] = (char *)calloc((end - start + 2),sizeof(char)))) { /*the last is 0 just to jump in some loops*/
					fprintf(file_output,"\nError: memory not reallocated. use_gff.23c \n");
					return 0; /*error*/
				}
				while(j<m && strcmp(fieldsgff[j].feature,"CDS")!=0) {j++;}
				seqid = fieldsgff[j].transcript_id;
				while(j<m && strcmp(fieldsgff[j].transcript_id,seqid) == 0) {
					ibeg = fieldsgff[j].start - start;
					iend = fieldsgff[j].end   - start;
					/*first check if CDS in the same transcript are overlapped: if yes discard*/
					l = ibeg;
					while(l<=iend && matrix_coding[i][l]==0) {l++;}
					if(l<=iend) {
						/*printf("\n Reading GTF file: CDS exons overlapped for the same transcript (gene_id = %s, transcript_id = %s). TRANSCRIPT NOT CONSIDERED.",fieldsgff[j].gene_id,fieldsgff[j].transcript_id);*/
						memset(matrix_coding[i],'\0',(end - start + 2)); /*reset the transcript to 0*/
						j=n;
						while(j<m) {
							while(strcmp(fieldsgff[j].feature,"CDS") != 0) {j++;}
							if(strcmp(fieldsgff[j].transcript_id,seqid) == 0) {
								strcpy(fieldsgff[j].feature,"CDS_EXCLUDED\0");
							} else break;
							j++;
						}
						break;
					}
					/*assign for each transcript the coding frame value to the matrix*/
					cstrand[0] = fieldsgff[j].strand[0];
					let[0] = fieldsgff[j].frame[0];
					if(cstrand[0]=='+') {
						for(l=ibeg;l<=iend;l++) {
							matrix_coding[i][l] = let[0];
							if(let[0]=='2') let[0] = '0';/*from 50 to 48*/
							else let[0]++;
						}
					}
					if(cstrand[0]=='-') {
						for(l=iend;l>=ibeg;l--) {
							matrix_coding[i][l] = let[0];
							if(let[0]=='2') let[0] = '0';/*from 48 to 50*/
							else let[0]++;
						}
					}
					j++;
					while(j<m && strcmp(fieldsgff[j].feature,"CDS")!=0) {j++;}
				}
			}
			/*discard the gene_id if the reading frames are not equal*/
			j=n;
			overlap_nf = 0;
			for(l=0;l<=end-start;l++) {
				i=0;
				while(i < ntransc && matrix_coding[i][l] == 0) {i++;}
				if(i<ntransc) 
					let[0] = matrix_coding[i][l];
				else 
					continue;
				while(i < ntransc && (matrix_coding[i][l] == 0 || matrix_coding[i][l] == let[0])) {i++;}
				if(i<ntransc) {overlap_nf = 1; break;}
			}
			
			if(overlap_nf==1) {
				/*printf("\n Reading GTF file: CDS exons with different reading frames or sense strand (gene_id = %s). GENE NOT CONSIDERED.",fieldsgff[j].gene_id);*/
				seqid = fieldsgff[n].gene_id;
				j=n;
				while(j<m) {
					while(strcmp(fieldsgff[j].feature,"CDS") != 0) {j++;}
					strcpy(fieldsgff[j].feature,"CDS_EXCLUDED\0");
					j++;
				}
			}
			else {
				/*use the max and min criteria to define fieldsgff2*/
				/*first define the final transcript in the row 'ntransc'*/
				if(!(matrix_coding[ntransc] = (char *)calloc((end - start + 2),sizeof(char)))) {
					fprintf(file_output,"\nError: memory not reallocated. use_gff.23b \n");
					return 0; /*error*/
				}
				
				if(strcmp(criteria_transcripts,"max")==0) {
					for(l=0;l<=end-start;l++) {
						i=0;
						while(i < ntransc && matrix_coding[i][l] == 0) {i++;} /*only count those positions with non-coding*/
						if(i<ntransc) /*if all transcripts are non-coding i==ntransc, so take the contrary*/
							matrix_coding[ntransc][l] = matrix_coding[i][l];
					}
				}
				if(strcmp(criteria_transcripts,"min")==0) {
					for(l=0;l<=end-start;l++) {
						i=0;
						while(i < ntransc && matrix_coding[i][l] != 0) {i++;} /*only count those positions and transcripts where coding exists*/
						if(i==ntransc)  /*if all transcripts have the position accept*/
							matrix_coding[ntransc][l] = matrix_coding[0][l];
					}
				}
				if(strcmp(criteria_transcripts,"first")==0) {
					for(l=0;l<=end-start;l++) {
						matrix_coding[ntransc][l] = matrix_coding[0][l];
					}
				}
				
				if(strcmp(criteria_transcripts,"long")==0) { /*look for the longest transcript*/
					longtr = (long int *) calloc(ntransc, sizeof(long int));
					for(i=0;i<ntransc;i++) {
						for(l=0;l<=end-start;l++) {
							if(matrix_coding[i][l] != 0) longtr[i] += 1;
						}
					}
					lotr = 0;
					for(i=1;i<ntransc;i++) {
						if(longtr[lotr] < longtr[i]) {
							lotr = i;
						}
					}
					free(longtr);
					for(l=0;l<=end-start;l++) {
						matrix_coding[ntransc][l] = matrix_coding[lotr][l];
					}
				}
				
				/*name of the consensus transcript_id field in all geneid*/
				sprintf(transcript,"%s_criteria_%s",criteria_transcripts,fieldsgff[j].gene_id);
				
				/*include the transcript in fieldsgff2: count all exons.*/
				
				l = 0;
				while(l <= end-start) {
					while(matrix_coding[ntransc][l] == 0) {l++;}
					jj = l;
					while(matrix_coding[ntransc][jj] != 0) {jj++;}
					if(jj-1<=end-start) {
						memcpy(fieldsgff2+nrows2,fieldsgff+j,sizeof(struct valuesgff)*1);/*include row in fieldsgff2*/
						fieldsgff2[nrows2].start     = l  + start;
						fieldsgff2[nrows2].end       = jj-1 + start;
						fieldsgff2[nrows2].strand[0] = cstrand[0];
						if(fieldsgff2[nrows2].strand[0] == '+') fieldsgff2[nrows2].frame[0]  = matrix_coding[ntransc][l];
						if(fieldsgff2[nrows2].strand[0] == '-') fieldsgff2[nrows2].frame[0]  = matrix_coding[ntransc][jj-1];
						strcpy(fieldsgff2[nrows2].feature,"CDS\0");
						strcpy(fieldsgff2[nrows2].transcript_id,transcript);
						nrows2++;
						/*
						if(nrows2 >= nrows) {
							if(!(fieldsgff2 = (struct valuesgff *)realloc(fieldsgff2,(nrows2+1)*sizeof(struct valuesgff)))) {
								fprintf(file_output,"\nError: memory not reallocated. use_gff.3b \n");
								return 0; 
							}
						}
						*/
					}
					l = jj;
				}
			}
			/*define all the rest of rows that are not CDS in fieldsgff2*/
			j = n;
			while(j<m) {
				while(j<m && strcmp(fieldsgff[j].feature,"CDS")==0) {j++;}
				if(j<m && strcmp(fieldsgff[j].feature,"CDS")!=0) {
					memcpy(fieldsgff2+nrows2,fieldsgff+j,sizeof(struct valuesgff)*1);/*include row in fieldsgff2*/
					nrows2++;
					/*
					if(nrows2 >= nrows) {
						if(!(fieldsgff2 = (struct valuesgff *)realloc(fieldsgff2,(nrows2+1)*sizeof(struct valuesgff)))) {
							fprintf(file_output,"\nError: memory not reallocated. use_gff.3b \n");
							return 0; 
						}
					}
					*/
					j++;
				}
			}
			
			/*free matrix_coding*/
			for(i=0;i<(ntransc+1);i++) free(matrix_coding[i]);
			free(matrix_coding);
			
			/*look for next gene_id*/
			n = m;
		}
		nrows = nrows2/*-1*/;
		/*reallocate rows for fieldsgff2*/
		if(!(fieldsgff2 = (struct valuesgff *)realloc(fieldsgff2,nrows*sizeof(struct valuesgff)))) {
			fprintf(file_output,"\nError: memory not reallocated. use_gff.3b \n");
			return 0; /*error*/
		}
		
		/*create a vector for erasing overlapped regions in different coding frame */
		if(!(vector_erase_overlapped = (int *)calloc(n_site,sizeof(int)))) {
			fprintf(file_output,"\nError: memory not reallocated. use_gff.23 \n");
			return 0; /*error*/
		}
		/*sort the struct by start to be fast in the comparison*/
		qsort(fieldsgff2,nrows/*-1*/,sizeof(struct valuesgff),comp_start_id);
		/*check that CDS regions with different gene_id are not overlapped: We assume we have replaced all '.' by its value in strand and frame */
		if(strcmp(subset_positions,"synonymous") == 0 || strcmp(subset_positions,"nonsynonymous") == 0 || strcmp(subset_positions,"silent") == 0 || strcmp(subset_positions,"0-fold") == 0 || strcmp(subset_positions,"2-fold") == 0 || strcmp(subset_positions,"3-fold") == 0 || strcmp(subset_positions,"4-fold") == 0) {
			for(n=0;n<nrows/*-1*/;n++) /*(n=0;n<nrows;n++)*/ {
				for(m=n+1;m<nrows;m++)/*(m=0;m<nrows;m++)*/ {
					if(m==n) continue;
					overlap_nf = 0;
					if(strcmp(fieldsgff2[m].feature,"CDS") == 0 && strcmp(fieldsgff2[n].feature,"CDS") == 0) {
						if(fieldsgff2[m].start >= fieldsgff2[n].start && fieldsgff2[m].start <= fieldsgff2[n].end) { /*overlapped*/
							/*check reading frames in strand*/
							if(fieldsgff2[m].strand[0] != fieldsgff2[n].strand[0]) overlap_nf = 1;
							else {
								/*check frame: LOOK FOR THE FRAME IN THE SAME POSITION (48 is ascii='0')*/
								if(fieldsgff2[n].strand[0]=='+') {
									rest  = fmod((double)(fieldsgff2[m].frame[0]),3.0);
									rest2 = fmod((double)(fieldsgff2[m].start-fieldsgff2[n].start + fieldsgff2[n].frame[0]),3.0);
								}
								if(fieldsgff2[n].strand[0]=='-') {
									rest  = fmod((double)(fieldsgff2[m].end-fieldsgff2[m].start + fieldsgff2[m].frame[0]),3.0);
									rest2 = fmod((double)(fieldsgff2[n].end-fieldsgff2[m].start + fieldsgff2[n].frame[0]),3.0);
								}
								if((int)rest != (int)rest2) 
									overlap_nf = 1;
							}
						}
						else {
							if(fieldsgff2[m].end >= fieldsgff2[n].start && fieldsgff2[m].end <= fieldsgff2[n].end) { /*overlapped*/
								/*check reading frames in  strand*/
								if(fieldsgff2[m].strand[0] != fieldsgff2[n].strand[0]) overlap_nf = 1;
								else {
									/*check frame: LOOK FOR THE FRAME IN THE SAME POSITION (48 is ascii='0')*/
									if(fieldsgff2[n].strand[0]=='+') {
										rest  = fmod((double)(fieldsgff2[m].end-fieldsgff2[m].start + fieldsgff2[m].frame[0]),3.0);
										rest2 = fmod((double)(fieldsgff2[m].end-fieldsgff2[n].start + fieldsgff2[n].frame[0]),3.0);
									}
									if(fieldsgff2[n].strand[0]=='-') {
										rest  = fmod((double)(fieldsgff2[m].frame[0]),3.0);
										rest2 = fmod((double)(fieldsgff2[n].end-fieldsgff2[m].end + fieldsgff2[n].frame[0]),3.0);
									}
									if((int)rest != (int)rest2) 
										overlap_nf = 1;
								}
							}
						}
					}
					if(overlap_nf == 1) {
						ibeg = fieldsgff2[m].start >= fieldsgff2[n].start ? fieldsgff2[m].start:fieldsgff2[n].start;
						iend = fieldsgff2[m].end   <= fieldsgff2[n].end ?   fieldsgff2[m].end:fieldsgff2[n].end;
						
						if(iend == fieldsgff2[m].end) {j = n; k = m;}
						else {j = m; k = n;/*ibeg == fieldsgff2[m].start*/}
						
						/*calculate the reading frame from the new starting point*/
						if(fieldsgff2[j].strand[0]=='+' && fieldsgff2[k].strand[0]=='+') {/*cut n*/
							cframe[0] = fieldsgff2[j].frame[0];
							rest = fmod((double)((fieldsgff2[j].frame[0]-48) + (iend-ibeg+1)),3.0);
							fieldsgff2[j].frame[0] = 48 + rest;
						}
						if(fieldsgff2[j].strand[0]=='-' && fieldsgff2[k].strand[0]=='-') {/*cut m*/
							cframe[0] = fieldsgff2[k].frame[0];
							rest = fmod((double)((fieldsgff2[k].frame[0]-48) + (iend-ibeg+1)),3.0);
							fieldsgff2[k].frame[0] = 48 + rest;
						}
						if(fieldsgff2[j].strand[0]=='+' && fieldsgff2[k].strand[0]=='-') {/*cut both*/
							cframe[0] = fieldsgff2[j].frame[0];
							rest = fmod((double)((fieldsgff2[j].frame[0]-48) + (iend-ibeg+1)),3.0);
							fieldsgff2[j].frame[0] = 48 + rest;
							cframe[0] = fieldsgff2[k].frame[0];
							rest = fmod((double)((fieldsgff2[k].frame[0]-48) + (iend-ibeg+1)),3.0);
							fieldsgff2[k].frame[0] = 48 + rest;
						}
						/*if(fieldsgff2[j].strand[0]=='-' && fieldsgff2[k].strand[0]=='+') no cut*/
						
						
						fieldsgff2[m].start = fieldsgff2[m].start >= fieldsgff2[n].start ? iend+1:fieldsgff2[m].start;
						fieldsgff2[n].start = fieldsgff2[n].start >= fieldsgff2[m].start ? iend+1:fieldsgff2[n].start;
						fieldsgff2[m].end = fieldsgff2[m].end <= fieldsgff2[n].end ? ibeg-1:fieldsgff2[m].end;
						fieldsgff2[n].end = fieldsgff2[n].end <= fieldsgff2[m].end ? ibeg-1:fieldsgff2[n].end;
						
						for(ii=ibeg;ii<iend;ii++) vector_erase_overlapped[ii] = 1;
						if(type_output == 0 || type_output == 10) {
							printf("\n Reading GTF file: Overlapping CDS regions with different reading frame (gene_id = %s vs gene_id = %s). OVERLAPPED REGION (from %ld to %ld) NOT CONSIDERED",fieldsgff2[m].gene_id,fieldsgff2[n].gene_id,ibeg+1,iend+1);
						}
					}
				}
			}
		}
		
		free(fieldsgff);
		
		/*write a file_gff with the rows be take into account*/
		*transcript = '\0';
		strcat(transcript,"_criteria_");
		strcat(transcript,criteria_transcripts);
		strcat(transcript,".gff");
		
		*name_fileinputgff2 = '\0';
		strncat(name_fileinputgff2,name_fileinputgff,256);
		ff = strstr(name_fileinputgff2,".");
		ff1 = ff;
		while(ff1 != NULL) {/*search the last '.'*/
			if((ff1 = strstr(ff+1,".")) != NULL)
				ff = ff1;
		}
		if(ff) {
			*ff = '\0';
			strcat(ff,transcript);
		}
		else strncat(name_fileinputgff,transcript,256);
		
		if((file_gff2 = fopen(name_fileinputgff2,"w")) == 0) {
			printf("\n  It is not possible to create the file %s",name_fileinputgff);
			return 0; 
		}
		fprintf(file_gff2,"#seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattributes\n");
		for(n=0;n<nrows;n++) {
			fprintf(file_gff2,"%s\t",fieldsgff2[n].filename);
			fprintf(file_gff2,"%s\t",fieldsgff2[n].source);
			fprintf(file_gff2,"%s\t",fieldsgff2[n].feature);
			fprintf(file_gff2,"%ld\t",fieldsgff2[n].start+1);
			fprintf(file_gff2,"%ld\t",fieldsgff2[n].end+1);
			fprintf(file_gff2,"%s\t",fieldsgff2[n].score);
			fprintf(file_gff2,"%c\t",fieldsgff2[n].strand[0]);
			if(fieldsgff2[n].frame[0]=='1') fprintf(file_gff2,"2\t");
			else {
				if(fieldsgff2[n].frame[0]=='2') fprintf(file_gff2,"1\t");
				else {
					fprintf(file_gff2,"%c\t",fieldsgff2[n].frame[0]);
				}
			}
			fprintf(file_gff2,"gene_id %s;",fieldsgff2[n].gene_id);
			fprintf(file_gff2,"transcript_id %s;",fieldsgff2[n].transcript_id);
			fprintf(file_gff2,"\n");
		}
		fclose(file_gff2);
		
		/*fieldsgff2, nrows and vector_erase_overlapped are the necessary information in this first part*/
		
		
		/**********************************************  END TO CHECK  *******************************************************/
		/**********************************************  END TO CHECK  *******************************************************/
		/**********************************************  END TO CHECK  *******************************************************/
		/**********************************************  END TO CHECK  *******************************************************/
		/**********************************************  END TO CHECK  *******************************************************/
		/**********************************************  END TO CHECK  *******************************************************/
		/**********************************************  END TO CHECK  *******************************************************/
		/**********************************************  END TO CHECK  *******************************************************/
		/*
		 
		 SECOND PART: ONCE THE GFF/GTF ROWS ARE INCLUDED AND FILTERED ACCORDING A OVERLAPPED TRANSCRIPTS CRITERIA, DETECT POSITIONS AND VARIANTS: 
		 
		 Warning: codons that belong to several transcripts will be counted consecutively in the -c max option. The codons that be cut will be discarded except for the first consecutive exon.
		 */
		
		/*counting frames already in annotation*/
		if(!(cframe_pos = (char *)calloc(n_site,sizeof(char)))) {
			fprintf(file_output,"\nError: memory not reallocated. use_gff.343 \n");
			free(fieldsgff);
			return 0; /*error*/
		}

		/*take the selected regions*/
		for(n=0;n<nrows;n++) {
			/*take seqname/gene_id of each row and check if there are coincidences before, if yes, not count*/
			seqid = fieldsgff2[n].seqname;/*if the same name, the exons must not be cut (if the frame is lost is not working)*/
			if(seqid[0]=='\0') seqid = fieldsgff2[n].gene_id;
			/*fileid = fieldsgff2[n].filename;*/
			i=n-1;
			k=0;
			while(i>=0) {
				if(strcmp(seqid,fieldsgff2[i].seqname) == 0 || strcmp(seqid,fieldsgff2[i].gene_id) == 0) {
					k = 1;
					break;
				}
				i--;
			}
			if(k) continue;
			
			/*include in current matrix_sizepos all rows with the same seqname/gene_id and the indicated feature*/	
			/*start frame and endframe are the start and the end position of the coding region for the seqid gene (whatever strand)*/
			
			if(strcmp(subset_positions,"synonymous") == 0 || strcmp(subset_positions,"nonsynonymous") == 0 || strcmp(subset_positions,"silent") == 0 || strcmp(subset_positions,"0-fold") == 0 || strcmp(subset_positions,"2-fold") == 0 || strcmp(subset_positions,"3-fold") == 0 || strcmp(subset_positions,"4-fold") == 0) {
				cstrand[0] = fieldsgff2[n].strand[0];
				/*frame is only read for the first CDS in the gene (at -/+ strand)*/
				cframe[0] = cframe[1] = 'N';
				startframe = endframe = -1;
				
				for(ii=0;ii<n_site;ii++) cmat[ii] = (double)0;
				if(strcmp(subset_positions,"silent") == 0) for(ii=0;ii<n_site;ii++) cmatnc[ii] = (double)1;
				
				for(j=n;j<nrows;j++) {
					if(strcmp(fieldsgff2[j].feature,"CDS") == 0) {
						if(strcmp(seqid,fieldsgff2[j].seqname) == 0 || strcmp(seqid,fieldsgff2[j].gene_id) == 0) {
							for(ii=fieldsgff2[j].start;ii<=fieldsgff2[j].end;ii++) cmat[ii] = (double)1;
							
							/*cframe_pos is used to check the continuos frame versus annotated frame*/
							if(fieldsgff2[j].strand[0] == '+') 
								cframe_pos[fieldsgff2[j].start] = fieldsgff2[j].frame[0];/*0 means continue, >0 means check*/
							if(fieldsgff2[j].strand[0] == '-') 
								cframe_pos[fieldsgff2[j].end] = fieldsgff2[j].frame[0];/*0 means continue, >0 means check*/

							if(cstrand[0] == 32 || cstrand[0] == '.') cstrand[0] = fieldsgff2[j].strand[0];
							else {
								if((cstrand[0] == '+' && fieldsgff2[j].strand[0] == '-') || (cstrand[0] == '-' && fieldsgff2[j].strand[0] == '+')) {
									cstrand[0] = '*';
									break;
								}
							}
							if(startframe == -1) {/*init*/
								cframe[0] = fieldsgff2[j].frame[0];
								cframe[1] = fieldsgff2[j].frame[0];
								startframe = fieldsgff2[j].start;
								endframe = fieldsgff2[j].end;
							}
							if(startframe > fieldsgff2[j].start) {
								startframe = fieldsgff2[j].start;
								cframe[0] = fieldsgff2[j].frame[0];
							}
							if(endframe < fieldsgff2[j].end) {
								endframe = fieldsgff2[j].end;
								cframe[1] = fieldsgff2[j].frame[0];
							}
						}
					}
					if(strcmp(subset_positions,"silent") == 0) {
						if(strcmp(fieldsgff2[j].feature,"CDS") == 0) {
							if(strcmp(seqid,fieldsgff2[j].seqname) == 0 || strcmp(seqid,fieldsgff2[j].gene_id) == 0)
								for(ii=fieldsgff2[j].start;ii<=fieldsgff2[j].end;ii++) cmatnc[ii] = (double)0;
						}
					}
				}
				if(cstrand[0] == 32) cstrand[0] = '+';
				if(cframe[0] == 32) cframe[0] = '0';
				if(cframe[1] == 32) cframe[1] = '0';
			}
			else {
				if(strcmp(subset_positions,"noncoding") == 0) {
					for(ii=0;ii<n_site;ii++) cmat[ii] = (double)1;
					for(j=n;j<nrows;j++) {
						if(strcmp(fieldsgff2[j].feature,"CDS") == 0) {
							if(strcmp(seqid,fieldsgff2[j].seqname) == 0 || strcmp(seqid,fieldsgff2[j].gene_id) == 0)
								for(ii=fieldsgff2[j].start;ii<=fieldsgff2[j].end;ii++) cmat[ii] = (double)0;
						}
					}
				}
				else {
					if(strcmp(subset_positions,"coding") == 0) {
						for(ii=0;ii<n_site;ii++) cmat[ii] = (double)0;
						for(j=n;j<nrows;j++) {
							if(strcmp(fieldsgff2[j].feature,"CDS") == 0) {
								if(strcmp(seqid,fieldsgff2[j].seqname) == 0 || strcmp(seqid,fieldsgff2[j].gene_id) == 0)
									for(ii=fieldsgff2[j].start;ii<=fieldsgff2[j].end;ii++) cmat[ii] = (double)1;
							}
						}
					}
					else {
						for(ii=0;ii<n_site;ii++) cmat[ii] = (double)0;
						for(j=n;j<nrows;j++) {
							if(strcmp(fieldsgff2[j].feature,subset_positions) == 0) {
								if(strcmp(seqid,fieldsgff2[j].seqname) == 0 || strcmp(seqid,fieldsgff2[j].gene_id) == 0)
									for(ii=fieldsgff2[j].start;ii<=fieldsgff2[j].end;ii++) cmat[ii] = (double)1;
							}
						}
					}
				}
			}
			/*count positions in matrix_sizepos*/
			/*first look at subset_positions if syn/nsyn/silent*/
			if(strcmp(subset_positions,"synonymous") == 0 || strcmp(subset_positions,"nonsynonymous") == 0 || strcmp(subset_positions,"silent") == 0 || strcmp(subset_positions,"0-fold") == 0 || strcmp(subset_positions,"2-fold") == 0 || strcmp(subset_positions,"3-fold") == 0 || strcmp(subset_positions,"4-fold") == 0) {
				/*search syn/nsyn positions: REMEMBER, LOOK AT 'CDS' FEATURE IN GTFv2 FORMAT!*/
				/*look at strand and frame fields*/
				if(cstrand[0] == '*') {
					if(strcmp(fieldsgff2[n].feature,"CDS") == 0) {
						if(type_output == 0 || type_output == 10) {
							/*printf("\n Error in GFF-format: annotation in %s is not considered. ",seqid);*/
							if(file_output) fprintf(file_output,"\n Error in GFF-format: annotation in %s is not considered. ",seqid);
						}
					}
					continue;
				}
				else {
					/*check frame*/
					if(cstrand[0] == '-') {
						start = (long int)n_site-1;
						end = 0;
						k = -1;
						endtrp = startframe;
						startframe = endframe;
						cframe[0] = cframe[1];
					}
					else {
						start = 0;
						end = (long int)n_site-1;
						endtrp = endframe;
						k = 1;
					}
					/*start looking at coding region sequences*/
					for(ii=start;ii*k <= end*k;ii += k) {							
						while(ii*k <= end*k && cmat[ii] == (double)0) ii += k;
						if(ii*k > end*k) break;
						if(ii == startframe) {
							if(cframe[0] == '1') {
								cmat[ii] = (double)0;
								do{
									ii += k;
								}while(ii*k <= end*k && cmat[ii] == (double)0);
								cmat[ii] = (double)0;
								continue;
							}
							if(cframe[0] == '2') {
								cmat[ii] = (double)0;
								continue;
							}
						}
						
						/*checking annotation versus continuous frame*/
						countframe = 0;
						ii2 = ii;
						if(cframe_pos[ii2] != 0) {
							if(cframe_pos[ii2] != '0') { /*start position of the codon*/
								if(cframe_pos[ii2] == '1') {
									cmat[ii] = (double)0;
									do{
										ii += k;
									}while(ii*k <= end*k && cmat[ii] == (double)0);
									if(ii*k > end*k) break;
									cmat[ii] = (double)0;
									continue;
								}
								if(cframe_pos[ii2] == '2') {
									cmat[ii] = (double)0;
									continue;
								}
							}
						}
						do{
							ii2 += k;
						}while(ii2*k <= end*k && cmat[ii2] == (double)0);
						if(ii2*k > end*k) break;
						countframe += 1;
						if(cframe_pos[ii2] != 0) {
							if(cframe_pos[ii2] != '1') { /*start position of the codon*/
								if(cframe_pos[ii2] == '2') {
									cmat[ii] = (double)0;
									do{
										ii += k;
									}while(ii*k <= end*k && cmat[ii] == (double)0);
									if(ii*k > end*k) break;
									cmat[ii] = (double)0;
									continue;
								}
								if(cframe_pos[ii2] == '0') {
									cmat[ii] = (double)0;
									continue;
								}
							}
						}
						do{
							ii2 += k;
						}while(ii2*k <= end*k && cmat[ii2] == (double)0);
						if(ii2*k > end*k) break;
						countframe += 1;
						if(cframe_pos[ii2] != 0) {
							if(cframe_pos[ii2] != '2') { /*start position of the codon*/
								if(cframe_pos[ii2] == '0') {
									cmat[ii] = (double)0;
									do{
										ii += k;
									}while(ii*k <= end*k && cmat[ii] == (double)0);
									if(ii*k > end*k) break;
									cmat[ii] = (double)0;
									continue;
								}
								if(cframe_pos[ii2] == '1') {
									cmat[ii] = (double)0;
									continue;
								}
							}
						}
						/*end checking annotation versus continuos frame*/
						
						/*function to read matrix with 3 * n_samp char, reverse-complementary if cstrand is '-', return 0 if more than 1 mutation in triplet. if gaps/uncertainty in one or two positions eliminate the triplet.*/
						if(tripletnsamp(cod3n,DNA_matr,cstrand[0],cmat,n_samp,(long int)n_site,end,ii,file_output/*,mainargc*/,include_unknown,type_output,nmhits,mhitbp,outgroup_presence,nsamoutg) == 0) {
							cmat[ii] = (double)0;
							do{
								ii += k;
							}while(ii*k <= end*k && cmat[ii] == (double)0);
							if(ii*k > end*k) break;
							cmat[ii] = (double)0;
							do{
								ii += k;
							}while(ii*k <= end*k && cmat[ii] == (double)0);
							if(ii*k > end*k) break;
							cmat[ii] = (double)0;
							continue;
						}
                        /************************** FOLD DEGENERATED POSITIONS: ONLY COUNT THE ANCESTRAL OR MORE FREQUENT *****/
                        if(strcmp(subset_positions,"0-fold") == 0 || strcmp(subset_positions,"2-fold") == 0 || strcmp(subset_positions,"3-fold") == 0 || strcmp(subset_positions,"4-fold") == 0) {
                            /*count biallelic sites in matrix_segrpos: REMEMBER overlapping coding regions are not well calculated!*/
                            jj = 0;hf1 = 1; hf2 = 0;jo=0;
                            while(cod3n[3*jj+0] == '5' || cod3n[3*jj+1] == '5'|| cod3n[3*jj+2] == '5') jj++;
                            for(j=jj+1;j<n_samp;j++) {
                                if(memcmp(cod3n+3*jj,cod3n+3*j,3) != 0) {
                                    if(cod3n[3*j+0] == '5' || cod3n[3*j+1] == '5'|| cod3n[3*j+2] == '5') continue;
                                    else {jo = j;hf2 += 1;}
                                } else hf1 += 1;
                            }
                            /*No ougroup: look for the higher frequency codon*/
                            if(hf1 < hf2) {
                                j = jo;
                            } else {
                                j = jj;
                                jj = jo;
                                
                            }
                            if(outgroup_presence) {
                                /*look for the outgroup codon*/
                                for(jo=n_samp-nsamoutg;jo<n_samp;jo++) {
                                    if(memcmp(cod3n+3*j,cod3n+3*jo,3) != 0) {
                                        if(cod3n[3*jo+0] == '5' || cod3n[3*jo+1] == '5'|| cod3n[3*jo+2] == '5') continue;
                                        else {
                                            jj = j;
                                            j = jo;
                                            break;}
                                    }
                                }
                            }
                            stop = 0;
                            cod3ft[0] = cod3ft[1] = cod3ft[2] = (double)0;
                            neffsam = 0;
                            for(q=0;q<64;q++) {
                                if(memcmp(tripletsN[q],cod3n+3*j,3) == 0) {
                                    aaseq[0] = genetic_code[q];
                                    break;
                                }
                            }
                            if(q==64 && include_unknown == 1)
                                continue;/*missing codon, not counting. It can be missing positions*/
                            if(aaseq[0] == '*') {
                                ii2 = ii;
                                do{
                                    ii2 += k;
                                }while(ii2*k <= end*k && cmat[ii2] == (double)0);
                                do{
                                    ii2 += k;
                                }while(ii2*k <= end*k && cmat[ii2] == (double)0);
                                if(endtrp != ii2) {
                                    if(type_output == 0 || type_output == 10) {
                                        /*if(mainargc > 1)*/ /*printf("\n Excluded codons:  Stop Codon starting at position %ld.",ii+1);*//**/
                                        if(file_output) fprintf(file_output,"\n Excluded codons: Stop Codon starting at position %ld.",ii+1);
                                    }
                                }
                                stop = 1;
                            }
                            countpath3 = 0;
                            for(i=0;i<3;i++) {
                                /*the number of positions at this codon for the selected positions*/
                                cod3f = (double)0;
                                countpath = 0;
                                if(strcmp(subset_positions,"0-fold") == 0) {
                                    countpath += 1;
                                    if(n_fold[q][i] == 0)
                                        cod3f += (double)1;
                                }
                                if(strcmp(subset_positions,"2-fold") == 0) {
                                    countpath += 1;
                                    if(n_fold[q][i] == 2)
                                        cod3f += (double)1;
                                }
                                if(strcmp(subset_positions,"3-fold") == 0) {
                                    countpath += 1;
                                    if(n_fold[q][i] == 3)
                                        cod3f += (double)1;
                                }
                                if(strcmp(subset_positions,"4-fold") == 0) {
                                    countpath += 1;
                                    if(n_fold[q][i] == 4)
                                        cod3f += (double)1;
                                }
                                if(countpath) {
                                    cod3ft[i] += cod3f/(double)countpath;
                                    countpath3 = 1;
                                }
                            }
                            if(countpath3)
                                neffsam += 1;
                        }
                        /**************************************************************************************************/
                        else {
                            /*count positions in cmat*/
                            stop = 0;
                            cod3ft[0] = cod3ft[1] = cod3ft[2] = (double)0;
                            neffsam = 0;
                            for(j=0;j<n_samp;j++) {
                                for(q=0;q<64;q++) {
                                    if(memcmp(tripletsN[q],cod3n+3*j,3) == 0) {
                                        aaseq[0] = genetic_code[q];
                                        break;
                                    }
                                }
                                if(q==64 && include_unknown == 1) 
                                    continue;/*missing codon, not counting. It can be missing positions*/
                                if(aaseq[0] == '*') {
                                    ii2 = ii;
                                    do{
                                        ii2 += k;
                                    }while(ii2*k <= end*k && cmat[ii2] == (double)0);
                                    do{
                                        ii2 += k;
                                    }while(ii2*k <= end*k && cmat[ii2] == (double)0);
                                    if(endtrp != ii2) {
                                        if(type_output == 0 || type_output == 10) {
                                            /*if(mainargc > 1)*/ /*printf("\n Excluded codons:  Stop Codon starting at position %ld.",ii+1);*//**/
                                            if(file_output) fprintf(file_output,"\n Excluded codons: Stop Codon starting at position %ld.",ii+1);
                                        }
                                    }
                                    stop = 1;
                                    break;
                                }
                                countpath3 = 0;
                                for(i=0;i<3;i++) {
                                    /*the number of positions at this codon for the selected positions*/
                                    cod3f = (double)0;
                                    countpath = 0;
                                    if((strcmp(subset_positions,"synonymous") == 0) ||
                                       (strcmp(subset_positions,"nonsynonymous") == 0)||
                                       (strcmp(subset_positions,"silent") == 0)) {
                                        for(nn=1;nn<=4;nn++) {
                                            /*do all possible combinations of one mutation in the codon to calculate the number of syn/nsyn positions*/
                                            cod3put[0] = cod3n[3*j+0];
                                            cod3put[1] = cod3n[3*j+1];
                                            cod3put[2] = cod3n[3*j+2];
                                            switch(nn) {
                                                case 1:
                                                    cod3put[i] = '1';
                                                    break;
                                                case 2:
                                                    cod3put[i] = '2';
                                                    break;
                                                case 3:
                                                    cod3put[i] = '3';
                                                    break;
                                                case 4:
                                                    cod3put[i] = '4';
                                                    break;
                                            }
                                            if(memcmp(cod3put,cod3n+3*j,3) == 0) continue;
                                            for(qq=0;qq<64;qq++) {
                                                if(memcmp(tripletsN[qq],cod3put,3) == 0) {
                                                    aaput[0] = genetic_code[qq];
                                                    break;
                                                }
                                            }
                                            
                                            if(aaput[0] == '*') continue;
                                            else countpath += 1;

                                            if(aaseq[0] == aaput[0]) {
                                                if((strcmp(subset_positions,"synonymous") == 0) ||
                                                   (strcmp(subset_positions,"silent") == 0)) {
                                                    cod3f += (double)1; /*syn*/
                                                }
                                            }else 
                                                if(strcmp(subset_positions,"nonsynonymous") == 0) {
                                                    cod3f += (double)1; /*nsyn*/
                                                }
                                        }
                                    }
                                    if(countpath) {
                                        cod3ft[i] += cod3f/(double)countpath;
                                        countpath3 = 1;
                                    }
                                }
                                if(countpath3)
                                    neffsam += 1;
                            }
						}
						if(stop == 1) {
							cmat[ii] = (double)0;
							do{
								ii += k;
							}while(cmat[ii] == (double)0);
							if(ii*k > end*k) break;
							cmat[ii] = (double)0;
							do{
								ii += k;
							}while(cmat[ii] == (double)0);
							if(ii*k > end*k) break;
							cmat[ii] = (double)0;
							continue;
						}
						/*count biallelic sites in matrix_segrpos: REMEMBER overlapping coding regions are not well calculated!*/
                        jj = 0;hf1 = 1; hf2 = 0;jo=0;
						while(cod3n[3*jj+0] == '5' || cod3n[3*jj+1] == '5'|| cod3n[3*jj+2] == '5') jj++;
						for(j=jj+1;j<n_samp;j++) {
							if(memcmp(cod3n+3*jj,cod3n+3*j,3) != 0) {
								if(cod3n[3*j+0] == '5' || cod3n[3*j+1] == '5'|| cod3n[3*j+2] == '5') continue;
                                else {jo = j;hf2 += 1;}
                            } else hf1 += 1;
						}
                        /*No ougroup: look for the higher frequency codon*/
                        if(hf1 < hf2) {
                            j = jo;
                        } else {
                            j = jj;
                            jj = jo;
                            
                        }
                        if(outgroup_presence) {
                            /*look for the outgroup codon*/
                            for(jo=n_samp-nsamoutg;jo<n_samp;jo++) {
                                if(memcmp(cod3n+3*j,cod3n+3*jo,3) != 0) {
                                    if(cod3n[3*jo+0] == '5' || cod3n[3*jo+1] == '5'|| cod3n[3*jo+2] == '5') continue;
                                    else {
                                        jj = j;
                                        j = jo;
                                        break;}
                                }
                            }
                        }
                        if(hf1 < n_samp-nsamoutg) {
							for(i=0;i<3;i++) if(memcmp(cod3n+3*jj+i,cod3n+3*j+i,1) != 0) break;/*Here (i) a difference!*/
                            for(qq=0;qq<64;qq++) {
								if(memcmp(tripletsN[qq],cod3n+3*jj,3) == 0) {
									aaseq[0] = genetic_code[qq];
									break;
								}
							}
							for(q=0;q<64;q++) {
								if(memcmp(tripletsN[q],cod3n+3*j,3) == 0) {
									aaput[0] = genetic_code[q];
									break;
								}
							}
							if(aaseq[0] == aaput[0]) {
								if((strcmp(subset_positions,"nonsynonymous") == 0)) {
									ii2 = ii;
									for(xx=0;xx<i;xx++) {
										do{
											ii2 += k;
										}while(ii2*k <= end*k && cmat[ii2] == (double)0); /*conditional: hf1 < n_samp-nsamoutg*/
									}
									matrix_segrpos[ii2] = (double)0; /*reject syn variant but keep position...*/
								}
							}
                            else {
                                if((strcmp(subset_positions,"synonymous") == 0) ||
                                   (strcmp(subset_positions,"silent") == 0)) {
                                    ii2 = ii;
                                    for(xx=0;xx<i;xx++) {
                                        do{
                                            ii2 += k;
                                        }while(ii2*k <= end*k && cmat[ii2] == (double)0);
                                    }
                                    matrix_segrpos[/*ii+i*k*/ii2] = (double)0; /*reject nsyn variant but keep position...*/
                                }
                            }
                            if((strcmp(subset_positions,"0-fold") == 0)) {
                                if(n_fold[q][i] != 0) { //ancestral or the higher frequency codon
                                    ii2 = ii;
                                    for(xx=0;xx<i;xx++) {
                                        do{
                                            ii2 += k;
                                        }while(ii2*k <= end*k && cmat[ii2] == (double)0);
                                    }
                                    matrix_segrpos[ii2] = (double)0; //reject non 0-fold variant but keep position...
                                }
                            }
                            if((strcmp(subset_positions,"2-fold") == 0)) {
                                if(n_fold[q][i] != 2) { //ancestral or the higher frequency codon
                                    ii2 = ii;
                                    for(xx=0;xx<i;xx++) {
                                        do{
                                            ii2 += k;
                                        }while(ii2*k <= end*k && cmat[ii2] == (double)0);
                                    }
                                    matrix_segrpos[ii2] = (double)0; //reject non 2-fold variant but keep position...
                                }
                            }
                            if((strcmp(subset_positions,"3-fold") == 0)) {
                                if(n_fold[q][i] != 3) { //ancestral or the higher frequency codon
                                    ii2 = ii;
                                    for(xx=0;xx<i;xx++) {
                                        do{
                                            ii2 += k;
                                        }while(ii2*k <= end*k && cmat[ii2] == (double)0);
                                    }
                                    matrix_segrpos[ii2] = (double)0; //reject non 3-fold variant but keep position...
                                }
                            }
                            if((strcmp(subset_positions,"4-fold") == 0)) {
                                if(n_fold[q][i] != 4) { //ancestral or the higher frequency codon
                                    ii2 = ii;
                                    for(xx=0;xx<i;xx++) {
                                        do{
                                            ii2 += k;
                                        }while(ii2*k <= end*k && cmat[ii2] == (double)0);
                                    }
                                    matrix_segrpos[ii2] = (double)0; //reject non 4-fold variant but keep position...
                                }
                            }
						}
						/*weight cmat positions*/
						cmat[ii] = cod3ft[0]/((double)neffsam/*n_samp*/);
						do{
							ii += k;
						}while(cmat[ii] == (double)0);
						cmat[ii] = cod3ft[1]/((double)neffsam/*n_samp*/);
						do{
							ii += k;
						}while(cmat[ii] == (double)0);
						cmat[ii] = cod3ft[2]/((double)neffsam/*n_samp*/);
					}
					/*include cmat in matrix_sizepos*/
					for(ii=0;ii<(long int)n_site;ii++) {
						if(matrix_sizepos[ii] < cmat[ii]) matrix_sizepos[ii] = cmat[ii];
					}
					/*noncoding for silent calculated separatedly*/
					if(strcmp(subset_positions,"silent") == 0) {
						for(ii=0;ii<(long int)n_site;ii++) {
							if(cmatsil[ii] > cmatnc[ii]) cmatsil[ii] = cmatnc[ii];
						}
					}
				}
			}
			else {
				/*include cmat in matrix_sizepos for other options*/
				if(strcmp(subset_positions,"noncoding") == 0) {
					for(ii=0;ii<(long int)n_site;ii++) {
						if(matrix_sizepos[ii] > cmat[ii]) matrix_sizepos[ii] = cmat[ii];
					}
				}
				else {
					for(ii=0;ii<(long int)n_site;ii++) {
						if(matrix_sizepos[ii] < cmat[ii]) matrix_sizepos[ii] = cmat[ii];
					}
				}
			}
		}
		/*in silent: add noncoding to syn positions*/
		if(strcmp(subset_positions,"silent") == 0) {
			for(ii=0;ii<(long int)n_site;ii++) {
				if(cmatsil[ii] == (double)1) matrix_sizepos[ii] = cmatsil[ii];
			}
		}
		
		/*erase overlapped regions with different reading frame*/
		for(ii=0;ii<(long int)n_site;ii++) {
			if(vector_erase_overlapped[ii] == 1) matrix_sizepos[ii] = 0;
		}
		
		free(cmat);
		free(cmatnc);
		free(cmatsil);
		if(strcmp(subset_positions,"synonymous") == 0 || strcmp(subset_positions,"nonsynonymous") == 0 || strcmp(subset_positions,"silent") == 0 || strcmp(subset_positions,"0-fold") == 0 || strcmp(subset_positions,"2-fold") == 0 || strcmp(subset_positions,"3-fold") == 0 || strcmp(subset_positions,"4-fold") == 0) {
			free(cod3n);
			free(cod3put);
		}
		free(fieldsgff2);
		free(vector_erase_overlapped);
		free(cframe_pos);
	}
	
	return 1; /*ok*/
}

/*function to read matrix with 3 * n_samp char, reverse-complementary if cstrand is '-', 
 if(include_unknown==0) return 0 if more than 1 mutation in triplet or gaps/uncertainty
 if(include_unknown==1) keep only the two more frequent variants in the codon, and allow missing*/
int tripletnsamp(char *cod3n,char *DNA_matr,char strand,double *cmat,
					int n_samp,long int n_site,long int end,long int ii,
					FILE *file_output/*,int mainargc*/,int include_unknown,int type_output,
					long int *nmhits, long int *mhitbp, int outgroup_presence, int nsamoutg)
{
	int i,j,k,jj,x,z;
	long int ii2;
	char triplet1[3],triplet2[3];
	char **codons_sam;
	int cbeg,cend;
	int *freq_codons;
	int variant1,variant2,fvar1,fvar2;
	int outgroup_codon = 0;
	
	if(strand == '-') k = -1;
	else k = 1;
    cbeg=cend=0;	
	/*make matrix with triplets*/
	ii2 = ii;
	for(i=0;i<3;i++)
	{
		for(j=0;j<n_samp;j++)
		{
			cod3n[3*j+i] = DNA_matr[n_site*j+ii2];
			if(DNA_matr[n_site*j+ii2] == '6') cod3n[3*j+i] = '5';
			/*read complementary when strand is '-'*/
			if(k==-1) {
				if(cod3n[3*j+i] == '1') cod3n[3*j+i] = '4';
				else {
					if(cod3n[3*j+i] == '2') cod3n[3*j+i] = '3';
					else {
						if(cod3n[3*j+i] == '3') cod3n[3*j+i] = '2';
						else if(cod3n[3*j+i] == '4') cod3n[3*j+i] = '1';
					}
				}
			}
		}

		if(i < 2) {
			do {
				ii2 += k;
			}while(ii2*k <= end*k && cmat[ii2] == (double)0);
		}

		if(ii2*k > end*k) 
			return 0;
	}
	
	if(include_unknown == 1) {
		/*
		 In case option -u=1, we "accept" positions with multiple mutations because we replace the lowest freqency (no outgroup) by 'N's
		 look for outgroup: if present, look for triplet in outgroup, if Ns, do it as no-outgroup
		 if outgroup and triplet, look for other triplets, keep only the more frequent triplet plus outgroup
		 if no outgroup keep the two more frequent triplet (modify DNA_matr)
		*/
		if((freq_codons = (int *)calloc(n_samp,sizeof(int))) == 0) {
			printf("\nError assigning memory: tripletnsamp.1a. Exit.\n");
			exit(0);
		} 
		if((codons_sam = (char **)calloc(n_samp,sizeof(char *))) == 0) {
			printf("\nError assigning memory: tripletnsamp.1b. Exit.\n");
			exit(0);
		} 
		for(j=0;j<n_samp;j++) {
			if((codons_sam[j] = (char *)calloc(3,sizeof(char))) == 0) {
				printf("\nError assigning memory: tripletnsamp.1c. Exit.\n");
				exit(0);
			} 
		}
		if(outgroup_presence == 1) {
			j = 0;
			outgroup_codon = 0;
			for(jj=n_samp-nsamoutg;jj<n_samp;jj++) { /*assign the first triplet from outgroup*/
				if(cod3n[3*jj+0] != '5' && cod3n[3*jj+1] != '5' && cod3n[3*jj+2] != '5') {
					codons_sam[j][0] = cod3n[3*jj+0];
					codons_sam[j][1] = cod3n[3*jj+1];
					codons_sam[j][2] = cod3n[3*jj+2];
					freq_codons[j] += 1;
					outgroup_codon = 1;
					break;
				}
			}
			if(codons_sam[j][0] == 0) { /*if the outgroup no exist*/
				for(jj=0;jj<n_samp-nsamoutg;jj++) { /*assign the first triplet from samples*/
					if(cod3n[3*jj+0] != '5' && cod3n[3*jj+1] != '5' && cod3n[3*jj+2] != '5') {
						codons_sam[j][0] = cod3n[3*jj+0];
						codons_sam[j][1] = cod3n[3*jj+1];
						codons_sam[j][2] = cod3n[3*jj+2];
						freq_codons[j] += 1;
						break;
					}
				}
			}
			if(codons_sam[j][0] == 0) {/*if the triplet no exist return zero values*/
				free(freq_codons);
				for(j=0;j<n_samp;j++) free(codons_sam[j]);
				free(codons_sam);		
				return 1; /*missing triplet*/
			}
			else j++;
		}
		else {
			j = 0;
			for(jj=0;jj<n_samp;jj++) { /*assign the first triplet from samples*/
				if(cod3n[3*jj+0] != '5' && cod3n[3*jj+1] != '5' && cod3n[3*jj+2] != '5') {
					codons_sam[j][0] = cod3n[3*jj+0];
					codons_sam[j][1] = cod3n[3*jj+1];
					codons_sam[j][2] = cod3n[3*jj+2];
					freq_codons[j] += 1;
					break;
				}
			}
			if(codons_sam[j][0] == 0) {/*if the triplet no exist return zero values*/
				return 1; /*missing triplet*/
			}
			else j++;
		}
		
		if(j==1) {
			/*calculate the frequencies of the triplet(s)*/
			cend = cbeg = 0;
			if(outgroup_presence == 1) {
				if(jj >= n_samp-nsamoutg) {cbeg = 0; cend = n_samp-nsamoutg;}
				if(jj <  n_samp-nsamoutg) {cbeg = jj+1; cend = n_samp-nsamoutg;}
			}
			else {cbeg = jj+1; cend = n_samp;}
			
			for(i=cbeg;i<cend;i++) {/*look at each sample*/
				if(cod3n[3*i+0] != '5' && cod3n[3*i+1] != '5' && cod3n[3*i+2] != '5') {/*look at valid codons*/
					for(k=0;k<j;k++) {/*look at each defined codon*/
						if(k==j-1 && memcmp(codons_sam[k],cod3n+(3*i),3*sizeof(char)) != 0) {
							codons_sam[j][0] = cod3n[3*i+0];
							codons_sam[j][1] = cod3n[3*i+1];
							codons_sam[j][2] = cod3n[3*i+2];
							freq_codons[j] += 1;
							j++;
							break;
						}
						else {
							if(memcmp(codons_sam[k],cod3n+(3*i),3*sizeof(char)) == 0) {
								freq_codons[k] += 1;
								break;
							}
						}
					}
				}
			}
			/*keep the highest frequency codon variants that be "compatible" (no more than one mutation per codon)*/
			/*the outgroup must be included in any case: outgroup_codon=1 if exist, otherwise 0*/
			/*keep no more than 2 variants, for the rest add misssing ('5')*/
			variant1 = variant2 = -1;
			fvar1 = fvar2 = 0;
			if(j>1) {
				if(j==2) {
					x=0; z=1;
					k = 0;
					if(codons_sam[x][0] == codons_sam[z][0]) k++;
					if(codons_sam[x][1] == codons_sam[z][1]) k++;
					if(codons_sam[x][2] == codons_sam[z][2]) k++;
					if(k<2) {/*multiple hit*/
						/*assign missing values to the other triplets*/
						for(i=0;i<n_samp;i++) {
							if(memcmp(codons_sam[z],DNA_matr+n_site*i+ii,3*sizeof(char)) == 0) {
								DNA_matr[n_site*i+ii+0] = '5';
								DNA_matr[n_site*i+ii+1] = '5';
								DNA_matr[n_site*i+ii+2] = '5';
							}
						}
						
						mhitbp[*nmhits] = ii;
						nmhits[0] = nmhits[0] + 1/*j-2*/;				
					}
				}
				else {
					if(outgroup_codon == 1) {
						variant1 = 0;
						fvar1 = freq_codons[0];
						variant2 = 1;
						fvar2 = freq_codons[1];
						for(k=2;k<j;k++) {
							if(fvar2 < freq_codons[k]) {
								variant2 = k;
								fvar2 = freq_codons[k];
							}
						}
					}
					else {
						variant1 = 0;
						fvar1 = freq_codons[0];
						for(k=1;k<j;k++) {
							if(fvar1 < freq_codons[k]) {
								variant1 = k;
								fvar1 = freq_codons[k];
							}
						}
						fvar2 = 0;
						for(k=0;k<j;k++) {
							if(k != variant1) {
								if(fvar2 < freq_codons[k]) {
									variant2 = k;
									fvar2 = freq_codons[k];
								}
							}
						}
					}
					/*assign missing values to the other triplets*/
					for(i=0;i<n_samp;i++) {
						for(k=0;k<j;k++) {
							if(k != variant1 && k != variant2) {
								if(memcmp(codons_sam[k],DNA_matr+n_site*i+ii,3*sizeof(char)) == 0) {
									DNA_matr[n_site*i+ii+0] = '5';
									DNA_matr[n_site*i+ii+1] = '5';
									DNA_matr[n_site*i+ii+2] = '5';
								}
							}
						}
					}
					mhitbp[*nmhits] = ii;
					nmhits[0] = nmhits[0] + 1/*j-2*/;				
				}
			}
			free(freq_codons);
			for(j=0;j<n_samp;j++) free(codons_sam[j]);
			free(codons_sam);
		}
	}
	if(include_unknown == 0) {
		/*check if gaps/uncertainty*/
		for(i=0;i<3;i++) {
			for(j=0;j<n_samp;j++) {
				if(!(cod3n[3*j+i] == '1' || cod3n[3*j+i] == '2' || cod3n[3*j+i] == '3' || cod3n[3*j+i] == '4' )) {
					/*if(mainargc > 1) printf(" GU %ld.",ii+1);*/
					if(cod3n[3*j+0] == '5' && cod3n[3*j+1] == '5' && cod3n[3*j+2] == '5') {
						if(type_output == 0 || type_output == 10) {
							/*printf("\n Excluded codons:  Gap/Uncertainty starting at position %ld.",ii+1);*/
							if(file_output) fprintf(file_output,"\n Excluded codons: Gap/Uncertainty starting at position %ld.",ii+1);
						}
					}
					return 0;
				}
			}
		}
		/*check if more than 1 mutation*/
		jj = 0;
		while(jj<n_samp && (cod3n[3*jj+0] == '5' || cod3n[3*jj+1] == '5'|| cod3n[3*jj+2] == '5')) jj++;
		triplet1[0] = cod3n[3*jj+0];
		triplet1[1] = cod3n[3*jj+1];
		triplet1[2] = cod3n[3*jj+2];
		triplet2[0] = '\0';
		for(i=jj+1;i<n_samp;i++) {
			if(memcmp(triplet1,cod3n+(3*i),3*sizeof(char)) != 0) { /*if different of triplet1*/
				if(triplet2[0] == '\0' || (cod3n[3*i+0] == '5' || cod3n[3*i+1] == '5'|| cod3n[3*i+2] == '5')) 
				{ /*triplet2 is undefined. we check for no missing values*/
					if(cod3n[3*i+0] != '5' && cod3n[3*i+1] != '5' && cod3n[3*i+2] != '5') {
						triplet2[0] = cod3n[3*i+0];
						triplet2[1] = cod3n[3*i+1];
						triplet2[2] = cod3n[3*i+2];					
						if(memcmp(triplet2,triplet1,3*sizeof(char)) != 0) {
							z=0; 
							for(x=0;x<3;x++) {
								if(triplet1[x]!=triplet2[x]) 
									z++;
							}
							if(z>1) { /*if more than one mutation --> mhit*/
								/*if(include_unknown == 0) {*/
									mhitbp[*nmhits] = ii;/*mhits are stablished by codons, not by positions*/
									nmhits[0] = nmhits[0] + 1;
									if(type_output == 0 || type_output == 10) {
										/*printf("\n Excluded codons: Multiple Mutations starting at position %ld.",ii+1);*/
										if(file_output) fprintf(file_output,"\n Excluded codons: Multiple Mutations starting at position %ld.",ii+1);
									}
									return 0;
								/*}*/
							}
						}
					}
				}
				else
				{	/*in case the triplet2 is defined, if we find a third triplet --> mhit*/
					if( memcmp(triplet2,cod3n+(3*i),3*sizeof(char)) != 0)
					{
						if(cod3n[3*i+0] != '5' && cod3n[3*i+1] != '5' && cod3n[3*i+2] != '5')
						{ /*check for no missing values*/
							/*if(include_unknown == 0) {*/
								mhitbp[*nmhits] = ii;/*mhits are stablished by codons, not by positions*/
								nmhits[0] = nmhits[0] + 1;
								if(type_output == 0 || type_output == 10) {
									/*if(mainargc > 1)*/ /*printf("\n Excluded codons: Multiple Mutations starting at position %ld.",ii+1);*/
									if(file_output) fprintf(file_output,"\n Excluded codons: Multiple Mutations starting at position %ld.",ii+1);
								}
								return 0;
							/*}*/
						}
					}
				}
			}
		}
	}
	
	return 1;
}

int comp_trcpt_id(const void *a,const void *b)
{
	int value;
	struct valuesgff *ia = (struct valuesgff *)a;
	struct valuesgff *ib = (struct valuesgff *)b;
	
	value = strcmp(ia->transcript_id,ib->transcript_id);
	return value;
}
int comp_start_id(const void *a,const void *b)
{
	struct valuesgff *ia = (struct valuesgff *)a;
	struct valuesgff *ib = (struct valuesgff *)b;
	
	if(ia->start < ib->start) return(-1);
	if(ia->start > ib->start) return(+1);
	return 0;
}
int comp_end_id(const void *a,const void *b)
{
	struct valuesgff *ia = (struct valuesgff *)a;
	struct valuesgff *ib = (struct valuesgff *)b;
	
	if(ia->end < ib->end) return(-1);
	if(ia->end > ib->end) return(+1);
	return 0;
}
int comp_gene_id(const void *a,const void *b)
{
	int value;
	struct valuesgff *ia = (struct valuesgff *)a;
	struct valuesgff *ib = (struct valuesgff *)b;
	
	value = strcmp(ia->gene_id,ib->gene_id);
	return value;
}

/*do a function that read the genetic code and add the degeneration to each position at each codon*/
int function_do_nfold_triplets(int n_fold[64][3], char *genetic_code, char tripletsN[64][3])
{
    int i,j,k,l;
    char check_tripletN[3];
    
    /*init*/
    for(i=0;i<64;i++) for(j=0;j<3;j++) n_fold[i][j] = 0;
    
    for(i=0;i<64;i++) {
        for(j=0;j<3;j++) {
            for(k='1';k<='4';k++) {
                if(tripletsN[i][j] != k) {
                    memcpy(check_tripletN, tripletsN[i],3*sizeof(char));
                    check_tripletN[j] = k;
                    l=0;
                    while(memcmp(check_tripletN, tripletsN[l], 3*sizeof(char)) && l<64) l++;
                    if(genetic_code[i] == genetic_code[l]) {
                        n_fold[i][j] += 1; /*degenerated*/
                    }
                    if(l==64)
                        return(0);/*error*/
                }
            }
        }
    }
    for(i=0;i<64;i++) for(j=0;j<3;j++) if(n_fold[i][j] > 0) n_fold[i][j] += 1;
    return(1);
}
