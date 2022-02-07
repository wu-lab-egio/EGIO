#include <stdio.h>
#include <string.h>

#define seqlen1 100000
#define seqlen2 100000
#define seqlen3 200000

char seq3_out[seqlen3] = "\0";
int scores[seqlen1][seqlen2] = { 0 };
int score_array[5][5] = { {5,-4,-4,-4,-2},{-4,5,-4,-4,-2},{-4,-4,5,-4,-2},{-4,-4,-4,5,-2},{-2,-2,-2,-2,-1} };

char* pairwisealign(char sequence1[], char sequence2[],char method[], int gap)
{
    void outputNW(char seq1[], char seq2[]);
    void scoresMatNW(char seq1[seqlen1], char seq2[seqlen2], char normol[], int gapopen);

    void outputSW(char seq1[], char seq2[]);
    void scoresMatSW(char seq1[seqlen1], char seq2[seqlen2], char normol[], int gapopen);

    char seq1[seqlen1] = "\0", seq2[seqlen2] = "\0", mtd[seqlen1] = "\0", normol[] = "ATCGN";
    strcpy(seq1, sequence1);
    strcpy(seq2, sequence2);
    strcpy(mtd, method);

    if (strcmp(mtd,"global")==0 | strcmp(mtd,"g")==0){
        scoresMatNW(seq1, seq2, normol, gap);
        outputNW(seq1, seq2);
    }
    else if(strcmp(mtd,"local")==0 | strcmp(mtd,"l")==0){
        scoresMatSW(seq1, seq2, normol, gap);
        outputSW(seq1, seq2);
    }
    else{
        printf("no method is called %s\n",mtd);
    }

    return seq3_out;
}

int max(int score1, int score2, int score3)
{
    if (score1 > score2)
        if (score1 > score3)
            return score1;
        else
            return score3;
    else
        if (score2 > score3)
            return score2;
        else
            return score3;
}

int find_index(char c, char normol[])
{
    int i = 0, index = 0;
    for (i = 0; i < strlen(normol); i++)
        if (c == normol[i])
            index = i;
    return index;
}

void scoresMatNW(char seq1[], char seq2[], char normol[], int gap)
{   
    int find_index(char c, char normol[]);
    int max(int score1, int score2, int score3);
    int index1 = 0, index2 = 0;
    int len1 = strlen(seq1), len2 = strlen(seq2);
    int score1 = 0, score2 = 0, score3 = 0, final_score = 0;
    int i = 0, j = 0, temp = 0;
    
    for (i = 1, temp = gap; i <= len1; i++, temp += gap)
        scores[i][0] = temp;
    for (j = 1, temp = gap; j <= len2; j++, temp += gap)
        scores[0][j] = temp;

    for(i=1;i<=len1;i++)
        for (j = 1; j <= len2; j++)
        {              
            index1 = find_index(seq1[i-1], normol);
            index2 = find_index(seq2[j-1], normol);
            score1 = scores[i-1][j] + gap;
            score2 = scores[i][j - 1] + gap;
            score3 = scores[i-1][j-1]+score_array[index1][index2];
            scores[i][j] = max(score1, score2, score3);
        }
}

void outputNW(char seq1[], char seq2[])
{   
    int max(int score1, int score2, int score3);
    char seq1_out[seqlen1] = "\0", seq2_out[seqlen2] = "\0";
    int score_max = 0, score1 = 0, score2 = 0, score3 = 0;
    int out_num = 0, seq1_num = strlen(seq1) - 1, seq2_num = strlen(seq2) - 1;
    int i, j, m;
    int acc = 0;
    int len1 = strlen(seq1), len2 = strlen(seq2);
    char seqstr[1] = {"|"};
    int start_i = 0, start_j = 0;
    int max_score = scores[len1][len2];

    // find the start position and fill gap

    if(len1 == len2)
    {
        start_i = strlen(seq1);
        start_j = strlen(seq2);
    }
    else if(len1 < len2)
    {
        for (j = len2; j > 0; j--)
            if (max_score <= scores[len1][j])
            {
                max_score = scores[len1][j];
                start_j = j;
            }
        start_i = len1;
        //---------
        if(start_j < len2)
        {
            for (j = len2; j > start_j; j--)
            {
                seq1_out[out_num] = '-';
                seq2_out[out_num] = seq2[j-1];
                out_num++;
            }
        }
    }
    else if(len1 > len2)
    {
        for (i = len1; i > 0; i--)
            if (max_score <= scores[i][len2])
            {
                max_score = scores[i][len2];
                start_i = i;
            }
        start_j = len2;
        //---------
        if(start_i < len1)
        {
            for (i = len1; i > start_i; i--)
            {
                seq1_out[out_num] = seq1[i-1];
                seq2_out[out_num] = '-';
                out_num++;
            }
        }
    }

    seq1_out[out_num] = seq1[start_i-1];
    seq2_out[out_num] = seq2[start_j-1];
    start_i = start_i - 1;
    start_j = start_j - 1;    
    out_num++;

    // backtrace: diagonal first
    while (start_i >= 1 && start_j >= 1)
    {   
        score1 = scores[start_i - 1][start_j];
        score2 = scores[start_i][start_j - 1];
        score3 = scores[start_i - 1][start_j - 1];
        score_max = max(score1, score2, score3);

        if (score3 == score_max)
        {
            seq1_out[out_num] = seq1[start_i-1];
            seq2_out[out_num] = seq2[start_j-1];
            start_i = start_i - 1;
            start_j = start_j - 1;
            out_num++;
        }
        else if (score1 == score_max)
        {
            
            seq1_out[out_num] = seq1[start_i-1];
            seq2_out[out_num] = '-';
            start_i = start_i - 1;
            out_num++;
        }
        else if (score2 == score_max)
        {
            
            seq2_out[out_num] = seq2[start_j-1];
            seq1_out[out_num] = '-';
            start_j = start_j - 1;
            out_num++;
        }
    }

    // fill the gap   
    if (start_i == 0 && start_j > 0)
    {
        while (start_j > 0)
        { 
            start_j = start_j - 1;
            seq2_out[out_num] = seq2[start_j];
            seq1_out[out_num] = '-';
            out_num++;
        }
    }
    else if (start_i > 0 && start_j == 0)
    {
        while (start_i > 0)
        {
            start_i = start_i - 1;
            seq1_out[out_num] = seq1[start_i];
            seq2_out[out_num] = '-';
            out_num++;
        }
    }

    //==========
    int out_len1 = strlen(seq1_out), out_len2 = strlen(seq2_out);
    for(i=out_len1-1;i>=0;i--)
    {
        seq3_out[acc]=seq1_out[i];
        acc ++;
    }
    for(m=0;m>=0;m--)
    {
        seq3_out[acc]=seqstr[m];
        acc++;
    }
    for(j=out_len2 - 1;j>=0;j--)
    {   
        seq3_out[acc]=seq2_out[j];
        acc++;
    }
    seq3_out[acc]='\0';
}

void scoresMatSW(char seq1[], char seq2[], char normol[], int gap)
{   
    int find_index(char c, char normol[]);
    int max(int score1, int score2, int score3);
    int index1 = 0, index2 = 0;
    int len1 = strlen(seq1), len2 = strlen(seq2);
    int score1 = 0, score2 = 0, score3 = 0, final_score = 0;
    int i = 0, j = 0, temp = 0;
    
    for (i = 1, temp = gap; i <= len1; i++, temp += gap)
        scores[i][0] = 0;
    for (j = 1, temp = gap; j <= len2; j++, temp += gap)
        scores[0][j] = 0;

    for(i=1;i<=len1;i++)
        for (j = 1; j <= len2; j++)
        {              
            index1 = find_index(seq1[i-1], normol);
            index2 = find_index(seq2[j-1], normol);
            score1 = scores[i-1][j] + gap;
            score2 = scores[i][j - 1] + gap;
            score3 = scores[i-1][j-1]+score_array[index1][index2];
            scores[i][j] = max(score1, score2, score3);

            if (scores[i][j] < 0)
            {
                scores[i][j] = 0;
            }
        }
}

void outputSW(char seq1[], char seq2[])
{   
    int max(int score1, int score2, int score3);
    char seq1_out[seqlen1] = "\0", seq2_out[seqlen2] = "\0";
    int score_max = 0, score1 = 0, score2 = 0, score3 = 0;
    int out_num = 0, seq1_num = strlen(seq1) - 1, seq2_num = strlen(seq2) - 1;
    int i, j, m;
    int acc = 0;
    int len1 = strlen(seq1), len2 = strlen(seq2);
    char seqstr[1] = {"|"};
    int start_i = strlen(seq1), start_j = strlen(seq2);
    int max_score = scores[0][0];

    for (i = 0; i <= len1; i++)
        for (j = 0; j <= len2; j++)
            if (max_score <= scores[i][j])
            {
                max_score = scores[i][j];
                start_j = j;
                start_i = i;
            }

    if (max_score > 0)
    {
        seq1_out[out_num] = seq1[start_i-1];
        seq2_out[out_num] = seq2[start_j-1];
        start_i = start_i - 1;
        start_j = start_j - 1;    
        out_num++;

        while (start_i >= 1 && start_j >= 1)
        {   
            score1 = scores[start_i - 1][start_j];
            score2 = scores[start_i][start_j - 1];
            score3 = scores[start_i - 1][start_j - 1];
            score_max = max(score1, score2, score3);

            if (score3 == score_max)
            {
                start_i = start_i - 1;
                start_j = start_j - 1;
                seq1_out[out_num] = seq1[start_i];
                seq2_out[out_num] = seq2[start_j];
                out_num++;
            }
            else if (score1 == score_max)
            {
                start_i = start_i - 1;
                seq1_out[out_num] = seq1[start_i];
                seq2_out[out_num] = '-';
                out_num++;
            }
            else if (score2 == score_max)
            {
                start_j = start_j - 1;
                seq2_out[out_num] = seq2[start_j];
                seq1_out[out_num] = '-';
                out_num++;
            }
        }
    }
    //==========
    int out_len1 = strlen(seq1_out), out_len2 = strlen(seq2_out);
    for(i=out_len1-1;i>=0;i--)
    {
        seq3_out[acc]=seq1_out[i];
        acc ++;
    }
    for(m=0;m>=0;m--)
    {
        seq3_out[acc]=seqstr[m];
        acc++;
    }
    for(j=out_len2 - 1;j>=0;j--)
    {   
        seq3_out[acc]=seq2_out[j];
        acc++;
    }
    seq3_out[acc]='\0';
}
