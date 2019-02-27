#include <stdio.h>
#include <stdlib.h>
#include <math.h>

unsigned int * XORSHIFT32(unsigned int seed, unsigned int n) // punctul 1
{
    unsigned int k, r;
    unsigned int *R = (unsigned int *)malloc(n*sizeof(unsigned int));
    r = seed;
    R[0] = seed;
    for(k = 1; k < n; k++)
    {
        r = r ^ (r << 13);
        r = r ^ (r >> 17);
        r = r ^ (r << 5);
        R[k] = r;
    }
    return R;
}

//structura pentru memorarea pixelilor
typedef struct
{
    unsigned int R, G, B;
}PIXEL;

int *read_fileBMP()
{
    //citirea de la tastatura caracter cu caracter
    char *file_pathBMP = (char *)malloc(sizeof(char));
    int dim_name_file_path=0;
    char c;
    do
    {
        c = getchar();
        char *aux = (char *)realloc(file_pathBMP,(dim_name_file_path+1)*sizeof(char));
        file_pathBMP = aux;
        if(c != '\n')
            file_pathBMP[dim_name_file_path++] = c;
    }while(c!='\n');

    file_pathBMP[dim_name_file_path] = NULL;

    //verificare pentru existenta fisierului
    FILE *f = fopen(file_pathBMP, "rb");
    if(f == NULL)
    {
        printf("Fisierul introdus nu exista!\n");
        return NULL;
    }

    fclose(f);
    return file_pathBMP;
}

PIXEL *linearizeArray_BMP(char *file_pathBMP) // punctul 2
{
    FILE *f = fopen(file_pathBMP, "rb");

    //dimensiunea fisierului (inaltime si latime)
    unsigned int widthW, heightH;
    fseek(f, 18, SEEK_SET);
    fread(&widthW, sizeof(unsigned int), 1, f);
    fread(&heightH, sizeof(unsigned int), 1, f);

    unsigned int dim = widthW * heightH ;
    PIXEL *c = (PIXEL *)malloc(dim * sizeof(PIXEL));

    //calcularea padding-ului
    int padding;
    if(widthW % 4 != 0)
        padding = 4 - (3 * widthW) % 4;
    else
        padding = 0;

    //liniarizarea pixelilor
    int i, j;
    unsigned char RGB[3];
    fseek(f, 54, SEEK_SET);
    for(i = heightH - 1; i >= 0; i--)
       {
           for(j = 0; j < widthW; j++)
            {
                fread(RGB, 3, 1, f);
                c[i * widthW + j].R = RGB[2];
                c[i * widthW + j].G = RGB[1];
                c[i * widthW + j].B = RGB[0];
            }
            fseek(f, padding, SEEK_CUR);
        }
    fclose(f);
    return c;
}

char * writeFile_linearizedArray(char *file_pathBMP, PIXEL *c) // punctul 3
{
    char *new_file = (char *)malloc(26 * sizeof(char));
    static unsigned int number = 0; //utilizat pentru a crea fisiere cu nume diferite la apelari diferite ale functiei
    if(number == 0) //pentru criptare
        {
            strcpy(new_file, "crypted_bmp.bmp");
            number = 1;
        }
    else
        if(number == 1) //pentru decriptare
        {
            strcpy(new_file, "decrypted_bmp.bmp");
            number = 2;
        }
        else
            if(number == 2) //pentru desenarea ferestrelor
            {
                strcpy(new_file, "image_with_pattern_detections.bmp");
            }
    FILE *fout = fopen(new_file, "wb");
    FILE *fin = fopen(file_pathBMP, "rb");

    //dimensiunea fisierului (lungime si latime)
    unsigned int widthW, heightH;
    fseek(fin, 18, SEEK_SET);
    fread(&widthW, sizeof(unsigned int), 1, fin);
    fread(&heightH, sizeof(unsigned int), 1, fin);

    //calcularea padding-ului imaginii
    int padding;
    if(widthW % 4 != 0)
        padding = 4 - (3 * widthW) % 4;
    else
        padding = 0;

    unsigned char header_byte;
    int i,j;

    rewind(fin);

    //citirea si scrierea header-ului dintr-un fisier in alt fisier
    for(i = 0; i < 54; i++)
    {
        fread(&header_byte, 1, 1, fin);
        fwrite(&header_byte, 1, 1, fout);
    }

    //scrierea imaginii liniarizate in vectorul c impreuna cu padding-ul
    unsigned char zero = '0';
    for(i = heightH - 1; i >= 0; i--)
    {
        for(j = 0; j < widthW; j++)
        {
            fwrite(&c[i * widthW + j].B, 1, 1, fout);
            fwrite(&c[i * widthW + j].G, 1, 1, fout);
            fwrite(&c[i * widthW + j].R, 1, 1, fout);
        }
        if(padding != 0)
            fwrite(&zero, padding * sizeof(unsigned char),1 ,fout);
    }

    fclose(fout);
    fclose(fin);

    //intoarcerea caii fisierului creat
    return new_file;
}


PIXEL * pixels_permutation_from_liniarizedArray(PIXEL *c, unsigned int *p, unsigned int n, unsigned int **pos) // vectorul pos este inversa permutarii
{
    unsigned int i;
    PIXEL *w = (PIXEL *)malloc(n * sizeof(PIXEL));
    *pos = (unsigned int *)malloc(n * sizeof(unsigned int));
    for(i = 0; i < n; i++)
        {
            w[p[i]] = c[i]; //permutarea
            (*pos)[p[i]] = i; //inversa permutarii
        }
    return w;
}

PIXEL XOR_pixel_pixel(PIXEL x, PIXEL y) //  XOR pixel cu pixel
{
    PIXEL a;
    a.B = x.B ^ y.B;
    a.G = x.G ^ y.G;
    a.R = x.R ^ y.R;
    return a;
}

//structura pentru memoraria si accesarea octetilor mai usor
typedef union
{
    unsigned int x;
    unsigned char bytes[4];
}bytes;


PIXEL XOR_pixel_number(PIXEL x, unsigned int num) //XOR pixel cu numar
{
    PIXEL a;
    bytes t = {num};
    a.B = x.B ^ t.bytes[0];
    a.G = x.G ^ t.bytes[1];
    a.R = x.R ^ t.bytes[2];
    return a;

}

PIXEL *array_Encryption(char *file_pathBMP, char **enc_file_path, char *path_secret_key) // punctul 4
{
    FILE *f = fopen(file_pathBMP, "rb");
    FILE *key = fopen(path_secret_key, "r"); //cheia secreta

    //dimensiunea fisierului (lungime si latime)
    unsigned int widthW, heightH;
    fseek(f, 18, SEEK_SET);
    fread(&widthW, sizeof(unsigned int), 1, f);
    fread(&heightH, sizeof(unsigned int), 1, f);

    //citirea cheii secrete S
    unsigned int secret_key, SV;
    fscanf(key, "%u", &secret_key);
    fscanf(key, "%u", &SV);

    //liniarizarea imaginii
    PIXEL *c = linearizeArray_BMP(file_pathBMP);

    unsigned int n = widthW * heightH;

    // algoritmul de xorshift - generare
    unsigned int *R = (unsigned int *)malloc((2 * n) * sizeof(unsigned int));
    unsigned int seed = secret_key;
    int i;
    R = XORSHIFT32(seed, 2 * n);

    //scrierea lui R intr-un fisier pentru utilizarea ulterioara in functia de decriptare
    FILE *g = fopen("file_R.txt", "wb");
    for(i = 0; i < 2 * n; i++)
        fwrite(&R[i], sizeof(unsigned int),1 , g);
    fclose(g);

    //algoritmul lui durstenfeld - permutare aleatorie
    unsigned int *p = (unsigned int *)malloc(n * sizeof(unsigned int));

    for(i = 0; i < n; i++)
        p[i] = i;

    for(i = n-1; i >= 1; i--)
    {
        unsigned int r = R[n-i] % (i+1);
        unsigned int aux = p[r];
        p[r] = p[i];
        p[i] = aux;
    }

    //permutarea pixelilor
    unsigned int *v;
    c = pixels_permutation_from_liniarizedArray(c, p, n, &v);

    //scrierea permutarii inverse intr-un fisier pentru utilizarea ulterioara in decriptare
    g = fopen("file_V.txt", "wb");
    for(i = 0; i < n; i++)
        fwrite(&v[i], sizeof(unsigned int),1 , g);
    fclose(g);

    PIXEL *ciphered_pixels = (PIXEL *)malloc(n * sizeof(PIXEL));

    //aplicarea criptarii conform cerintei
    PIXEL a;
    a = XOR_pixel_number(c[0], R[n]);
    a = XOR_pixel_number(a, SV);
    ciphered_pixels[0] = a;

    for(i = 1; i < n; i++)
        {
            a = XOR_pixel_pixel(c[i], ciphered_pixels[i - 1]);
            a = XOR_pixel_number(a, R[n + i]);
            ciphered_pixels[i] = a;
        }

    //scrierea in fisier a imaginii criptate si salvarea caii imaginii criptate
    *enc_file_path = writeFile_linearizedArray(file_pathBMP, ciphered_pixels);

    //eliberarea memoriei
    free(v);
    free(c);
    free(R);

    return ciphered_pixels;
}


PIXEL * array_Decryption(char *file_pathBMP, char *enc_file_path, char *path_secret_key) // punctul 5
{
    FILE *f = fopen("file_R.txt", "rb");
    FILE *g = fopen("file_V.txt", "rb");
    FILE *h = fopen(file_pathBMP, "rb");
    FILE *key = fopen(path_secret_key, "r");
    rewind(f); rewind(g); rewind(key);

    //dimensiunea fisierului (lungime si latime)
    unsigned int widthW, heightH;
    fseek(h, 18, SEEK_SET);
    fread(&widthW, sizeof(unsigned int), 1, h);
    fread(&heightH, sizeof(unsigned int), 1, h);
    fclose(h);

    //citirea cheii secrete
    unsigned int secret_key, SV;
    fscanf(key, "%u", &secret_key);
    fscanf(key, "%u", &SV);
    fclose(key);

    unsigned int n = widthW * heightH;
    int i;

    //citirea lui r din fisier si stergerea fisierului din memorie
    unsigned int *R = (unsigned int *)malloc((2 * n) * sizeof(unsigned int));
    for(i = 0; i < 2*n; i++)
        fread(&R[i], sizeof(unsigned int), 1, f);
    fclose(f); remove("file_R.txt");

    //citirea permutarii inverse din fisier si stergerea fisierului din memorie
    unsigned int *p = (unsigned int *)malloc(n * sizeof(unsigned int));
    for(i = 0; i < n; i++)
        fread(&p[i], sizeof(unsigned int), 1, g);
    fclose(g); remove("file_V.txt");

    //citirea pixelilor criptati cu reutilizare a functiei liniarize Array_BMP
    PIXEL *ciphered_pixels = (PIXEL *)malloc(n * sizeof(PIXEL));
    ciphered_pixels = linearizeArray_BMP("crypted_bmp.bmp");

    //aplicarea decriptarii
    PIXEL *deciphered_pixels = (PIXEL *)malloc(n * sizeof(PIXEL));

    PIXEL a;
    a = XOR_pixel_number(ciphered_pixels[0], R[n]);
    a = XOR_pixel_number(a, SV);
    deciphered_pixels[0] = a;

    for(i = 1; i < n; i++)
    {
        a = XOR_pixel_pixel(ciphered_pixels[i], ciphered_pixels[i - 1]);
        a = XOR_pixel_number(a, R[n + i]);
        deciphered_pixels[i] = a;
    }
    unsigned int *pos;

    //aplicarea inversei permutarii
    deciphered_pixels = pixels_permutation_from_liniarizedArray(deciphered_pixels, p, n, &pos);

    //scrierea in fisier a imaginii
    *enc_file_path = writeFile_linearizedArray("crypted_bmp.bmp", deciphered_pixels);

    //eliberearea memoriei
    free(R);
    free(p);
    free(ciphered_pixels);
    free(pos);

    return deciphered_pixels;
}

void chi(char *file_pathBMP)
{
    FILE *f = fopen(file_pathBMP, "rb");

    //dimensiunea fisierului (lungime + latime)
    unsigned int widthW, heightH;
    fseek(f, 18, SEEK_SET);
    fread(&widthW, sizeof(unsigned int), 1, f);
    fread(&heightH, sizeof(unsigned int), 1, f);

    unsigned int n = widthW * heightH;

    //calcularea padding-ului
    int padding;
    if(widthW % 4 != 0)
        padding = 4 - (3 * widthW) % 4;
    else
        padding = 0;

    //alocarea memoriei dinamic cu valori de 0 pentru a calcula frecventa culorilor RGB
    unsigned long *frecR = (unsigned long *)calloc(256, sizeof(unsigned long));
    unsigned long *frecG = (unsigned long *)calloc(256, sizeof(unsigned long));
    unsigned long *frecB = (unsigned long *)calloc(256, sizeof(unsigned long));

    //calcularea valorilor pentru fiecare pixel de culoare
    unsigned char RGB[3];
    int i, j;
    fseek(f, 54, SEEK_SET);
    for(i = 0; i < heightH; i++)
    {
        for(j = 0; j < widthW; j++)
        {
            fread(RGB, 3, 1, f);
            frecR[RGB[2]]++;
            frecG[RGB[1]]++;
            frecB[RGB[0]]++;
        }
        fseek(f, padding, SEEK_CUR);
    }

    double f_bar = n / 256;

    //calcularea valorile dintre 0 si 255 pe fiecare culoare in parte
    double chiR, chiG, chiB;
    for(i = 0; i<= 255; i++)
    {
        chiR = chiR + (((frecR[i] - f_bar)*(frecR[i] - f_bar))/f_bar);
        chiG = chiG + (((frecG[i] - f_bar)*(frecG[i] - f_bar))/f_bar);
        chiB = chiB + (((frecB[i] - f_bar)*(frecB[i] - f_bar))/f_bar);
    }

    printf("(%.2Lf, %.2Lf, %.2Lf)\n", chiR, chiG, chiB);

    //eliberearea memoriei
    free(frecR);
    free(frecG);
    free(frecB);
}


//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//
//------------------------------------------------------------------------------//

//functii ce tin de partea a doua din enunt

void grayscale_image(char* nume_fisier_sursa,char* nume_fisier_destinatie) //functia ce transforma o imagine color in nuante de gri
{
   FILE *fin, *fout;
   unsigned int dim_img, latime_img, inaltime_img;
   unsigned char pRGB[3], header[54], aux;

   fin = fopen(nume_fisier_sursa, "rb");

   fout = fopen(nume_fisier_destinatie, "wb+");

   fseek(fin, 18, SEEK_SET);
   fread(&latime_img, sizeof(unsigned int), 1, fin);
   fread(&inaltime_img, sizeof(unsigned int), 1, fin);

   //copiaza octet cu octet imaginea initiala in cea noua
	fseek(fin,0,SEEK_SET);
	unsigned char c;
	while(fread(&c,1,1,fin)==1)
	{
		fwrite(&c,1,1,fout);
		fflush(fout);
	}
	fclose(fin);

	//calculam padding-ul pentru o linie
	int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;

	fseek(fout, 54, SEEK_SET);
	int i,j;
	for(i = 0; i < inaltime_img; i++)
	{
		for(j = 0; j < latime_img; j++)
		{
			//citesc culorile pixelului
			fread(pRGB, 3, 1, fout);
			//fac conversia in pixel gri
			aux = 0.299*pRGB[2] + 0.587*pRGB[1] + 0.114*pRGB[0];
			pRGB[0] = pRGB[1] = pRGB[2] = aux;
        	fseek(fout, -3, SEEK_CUR);
        	fwrite(pRGB, 3, 1, fout);
        	fflush(fout);
		}
		fseek(fout,padding,SEEK_CUR);
	}
	fclose(fout);
}

//structura pentru pozitia fiecarei ferestre (pox_x, pos_y), impreuna cu culorile RGB si corelatia ei
typedef struct
{
    unsigned int pos_x, pos_y, R, G, B;
    double corr;
}window_pos;

// media valorilor intensitatii grayscale
double pixels_intensity_S(PIXEL *c, unsigned int widthW, unsigned int heightH)
{
    int i, j;
    double S = 0;
    for(i = heightH - 1; i >= 0; i--)
        for(j = 0; j < widthW; j++)
            S = S + c[i * widthW + j].R;
    S = S / 165;
    return S;
}

// deviatia standard a valorilor intensitatilor grayscale a pixelilor in sablonul S
double deviation_S(PIXEL *c, unsigned int widthW, unsigned int heightH, double S_bar, unsigned int n)
{

    int i, j;
    double S, sum = 0;
    for(i = heightH - 1 ; i >= 0; i--)
        for(j = 0; j < widthW; j++)
        {
            S = c[i * widthW + j].R;
            sum = sum + (S - S_bar) * (S - S_bar);

        }
    sum = sqrt(sum / (n - 1));
    return sum;
}

// media valorilor intensitatii grayscale a pixelilor din fereastra f
double pixel_intensity_fi(PIXEL *c, unsigned int image_widthW, unsigned int number_widthW, unsigned int number_heightH, int x, int y)
{
    int i, j;
    double s = 0;
    for(i = x + number_heightH - 1; i >= x; i--)
        for(j = y; j < y + number_widthW; j++)
            s = s + c[i * image_widthW + j].R;

    s = s / 165;
    return s;
}

// deviatia standard a valorilor intensitatilor grayscale a pixelilor in fereastra f
double deviation_fi(PIXEL *c, unsigned int image_widthW, unsigned int number_widthW, unsigned int number_heightH, double fi_bar, int x, int y, int n)
{
    int i, j;
    double fi, sum = 0;
    for(i = x + number_heightH - 1; i >= x; i--)
        for(j = y; j < y + number_widthW; j++)
        {
            fi = c[i * image_widthW + j].R;
            sum = sum + (fi - fi_bar) * (fi - fi_bar);
        }
    sum = sqrt( sum / (n - 1));
    return sum;
}

window_pos * template_match(PIXEL *c, char *image_path, char *number_template, float Ps, unsigned int *size_v, unsigned int *size_d, window_pos **D, unsigned int colorR, unsigned int colorG, unsigned int colorB) // punctul 7
{
    FILE *f = fopen(image_path, "rb"); //imagine
    FILE *g = fopen(number_template, "rb"); //sablon

    //dimensiunile sablonului (lungime + latime)
    unsigned int number_widthW, number_heightH;
    fseek(g, 18, SEEK_SET);
    fread(&number_widthW, sizeof(unsigned int), 1, g);
    fread(&number_heightH, sizeof(unsigned int), 1, g);
    fclose(g);

    //dimensiunile imaginii (lungime + latime)
    unsigned int image_widthW, image_heightH;
    fseek(f, 18, SEEK_SET);
    fread(&image_widthW, sizeof(unsigned int), 1, f);
    fread(&image_heightH, sizeof(unsigned int), 1, f);
    fclose(f);

    unsigned int n = number_widthW * number_heightH;

    //vectorul v este folosit pentru memorarea ferestrelor
    window_pos *v = (window_pos *)malloc(sizeof(window_pos));
    unsigned int size = 0;
    int i, j, k, l;

    //liniarizarea pixelilor din sablon
    PIXEL *h = linearizeArray_BMP(number_template);

    double  S, S_bar, fi, fi_bar, dev_s, dev_fi;

    S_bar = pixels_intensity_S(h, number_widthW, number_heightH);
    dev_s = deviation_S(h, number_widthW, number_heightH, S_bar, n);

    for(i = 0; i < image_heightH && i + number_heightH < image_heightH; i++) //cele doua for-uri cu i si j "merg" prin imagine
        for(j = 0; j < image_widthW && j + number_widthW < image_widthW; j++)
        {
            S = 0;
            fi_bar = pixel_intensity_fi(c, image_widthW, number_widthW, number_heightH, i, j);
            dev_fi = deviation_fi(c, image_widthW, number_widthW, number_heightH, fi_bar, i, j, n);

            for(k = i; k < i + number_heightH; k++) //cele doua for-uri cu k si l "merg" prin sablon
                for(l = j; l < j + number_widthW; l++)
                    S = S + (c[k * image_widthW + l].R - fi_bar) * ( h[(k - i) * number_widthW + (l - j)].R - S_bar ); //calcularea sumei din formula corelatiei

            //calcularea corelatiei conform formulei
            S = (double) S  / ( n * dev_s * dev_fi);

            if(S >= Ps) // verificarea corelatiei daca ea este mai mare decat pragul impus Ps = 0.5
            {

                //alocarea spatiului pentru inca o fereastra
                window_pos *aux = (window_pos *)realloc(v, (size + 1) * sizeof(window_pos));
                v = aux;

                //memorarea valorilor ferestrei curente
                v[size].pos_x = i;
                v[size].pos_y = j;
                v[size].R = colorR;
                v[size].G = colorG;
                v[size].B = colorB;
                v[size++].corr = S;


                //memorarea detectiilor, pe rand, pentru toate cele 10 sabloane
                window_pos *aux2 = (window_pos *)realloc(*D, (*size_d + 1) * sizeof(window_pos));
                *D = aux2;
                (*D)[*size_d] = v[size-1];
                (*size_d)++;
            }
        }

    //dimensiunile vectorilor v si D
    *size_v = size - 1;
    *size_d = *size_d - 1;
    return v;
}

//colorarea unei ferestre
PIXEL * colorize_window(PIXEL *c, int x, int y, unsigned int number_widthW, unsigned int number_heightH, unsigned int image_widthW, PIXEL color) //punctul 8
{
    int i, j;
    j = y;
    //modificarea pixelilor de pe cele doua laturi care corespund inaltimii, conform schemei de culoare impuse
    for(i = x; i < x + number_heightH; i++)
    {
        c[i * image_widthW + j] = color; //latura 1 (cea din stanga)
        c[i * image_widthW + (j + number_widthW -1)] = color; //latura 2 (cea din dreapta)
    }

     i = x;

     //modificarea pixelilor de pe cele doua laturi care corespund latimii, conform schemei de culoare impuse
     for(j = y; j < y + number_widthW; j++)
     {
        c[i * image_widthW + j] = color; //latura 1 (cea de sus)
        c[(i + number_heightH - 1) * image_widthW + j] = color; //latura 2 (cea de jos)
     }
     return c;
}

//colorarea tuturor ferestrelor
PIXEL * colorize_allWindows(PIXEL *c, char *image_path, char *number_template, window_pos *v, unsigned int size)
{
    FILE *f = fopen(image_path, "rb"); //pentru imagine
    FILE *g = fopen(number_template, "rb"); //pentru sablon

    //dimensiunea sablonului (inaltime + latime)
    unsigned int number_widthW, number_heightH;
    fseek(g, 18, SEEK_SET);
    fread(&number_widthW, sizeof(unsigned int), 1, g);
    fread(&number_heightH, sizeof(unsigned int), 1, g);
    fclose(g);

    //dimensiunea imaginii (inaltime + latime)
    unsigned int image_widthW, image_heightH;
    fseek(f, 18, SEEK_SET);
    fread(&image_widthW, sizeof(unsigned int), 1, f);
    fread(&image_heightH, sizeof(unsigned int), 1, f);
    fclose(f);

    int i, j;
    //colorarea fiecarei ferestre conform schemei de culoare
    for(i = 0; i < size; i++)
        {
            PIXEL color;
            color.R = v[i].R;
            color.G = v[i].G;
            color.B = v[i].B;
            c = colorize_window(c, v[i].pos_x, v[i].pos_y, number_widthW, number_heightH, image_widthW, color);
        }

    return c;
}

//functia de comparare pentru qsort
int cmp(const void *a, const void *b) // punctul 9
{
    window_pos va = *(window_pos *)a;
    window_pos vb = *(window_pos *)b;
    if(va.corr > vb.corr)
        return -1;
    else
        if(va.corr < vb.corr)
            return 1;
        else
            return 0;
}

//functia pentru verificarea a doua ferestre daca se suprapun si corespund parametrilor ceruti
int window_overlay(int x1, int y1, int x2, int y2, unsigned int number_widthW, unsigned int number_heightH)
{
    double s;
    unsigned int l, h, A, intersection;
    A = number_heightH * number_widthW; //aria sablonului

    //cazuri suprapunere
    if(x1 < x2 && y1 < y2) //caz 1
    {
        if((x2 - x1 >= number_heightH) || (y2 - y1 >= number_widthW)) //cand nu se suprapune
            return 0;

        l = y1 + number_widthW - y2;
        h = x1 + number_heightH - x2 ;
    }
    else
        if(x1 < x2 && y1 == y2) //caz 2
        {
            if(x2 - x1 >= number_heightH)
                return 0;

            l = number_widthW;
            h = x1 + number_heightH - x2;

        }
        else
            if(x1 < x2 && y1 > y2) //caz 3
            {
                if((x2 - x1 >= number_heightH) || (y1 - y2 >= number_widthW))
                    return 0;
                l = y2 + number_widthW - y1;
                h = x1 + number_heightH - x2;
            }
            else
                if(x1 == x2 && y1 > y2) //caz 4
                {
                    if(y1 - y2 >=number_widthW)
                        return 0;
                    l = y2 + number_widthW - y1;
                    h = number_heightH;
                }
                else
                    if(x1 > x2 && y1 > y2) //caz 5
                    {
                        if((x1 - x2 >= number_heightH) || (y1 - y2 >= number_widthW))
                            return 0;
                        l = y2 + number_widthW - y1;
                        h = x2 + number_heightH - x1;
                    }
                    else
                        if(x1 > x2 && y1 == y2) //caz 6
                        {
                            if(x1 - x2 >= number_heightH)
                                return 0;
                            l = number_widthW;
                            h = x2 + number_heightH - x1;
                        }
                        else
                            if(x1 > x2 && y1 < y2) //caz 7
                            {
                                if((x1 - x2 >= number_heightH) || (y2 - y1 >= number_widthW))
                                    return 0;
                                l = y1 + number_widthW - y2;
                                h = x2 + number_heightH - x1;
                            }
                            else
                                if(x1 == x2 && y1 < y2) //caz 8
                                {
                                    if(y2 - y1 >= number_widthW)
                                        return 0;
                                    l = y1 + number_heightH - y2;
                                    h = number_heightH;
                                }
                                else
                                    if(x1 == x2 && y1 == y2) //caz 9
                                    {
                                        l = number_widthW;
                                        h = number_widthW;
                                    }

    intersection = l * h; //aria intersectiei

    s = (double)(intersection) / (2 * A - intersection); //suprapunerea
    if(s > 0.2)
        return 1;
    return 0;

}

//stergerea unei valori care are corelatia necorespunzatoare
void delete_oneValue(window_pos **D, unsigned int *num, int pos)
{
    int i;
    for(i = pos; i < (*num) - 1; i++)
        (*D)[i] = (*D)[i + 1];
    (*num)--;
    *D = (window_pos *)realloc(*D, (*num) * sizeof(window_pos));

}

//stergerea tuturor valorile care au corelatiile necorespunzatoare cu ferestrele care se suprapun
window_pos * delete_allValues(window_pos *D, unsigned int number_widthW, unsigned int number_heightH, unsigned int *num)
{
    int i, j;
    for(i = 0; i < *num - 1; i++) //parcurgerea vectorului D si verificarea pentru fiecare pereche (i,j) cu i<j
        {
            for(j = i + 1; j < *num; j++)
                if(window_overlay(D[i].pos_x, D[i].pos_y, D[j].pos_x, D[j].pos_y, number_widthW, number_heightH) == 1) //daca se suprapun
                     if(D[i].corr > D[j].corr) //daca are o corelatie mai buna
                        {
                            unsigned int nr = *num;
                            delete_oneValue(&D, &nr, j);
                            *num = nr;
                            j--;
                        }
        }
    return D;
}

int main()
{
    //prima parte din enunt
    printf("Introdu calea fisierului ce urmeaza criptata: ");
    char *file_pathBMP = read_fileBMP(); //citirea de la tastatura a caii imaginii

    if(file_pathBMP != NULL) //sa existe calea fisierului introdus
    {
        printf("\nIntrodu calea cheii secrete: ");
        char *path_secret_key = read_fileBMP(); //citirea de la tastatura a caii cheii secrete

        if(path_secret_key != NULL) // sa existe calea pentru cheia secreta
        {
            PIXEL *c, *v;
            char *enc_file_path;

            //criptarea imaginii
            array_Encryption(file_pathBMP, &enc_file_path, path_secret_key);

            //decriptarea imaginii
            array_Decryption(file_pathBMP, enc_file_path, path_secret_key);

            printf("\nTestul chi pe imaginea initiala:\n");
            chi(file_pathBMP);

            printf("\nTestul chi pe imaginea criptata: \n");
            chi("crypted_bmp.bmp");
            printf("\n");
        }
        else //cazul de inexistenta
        {
            exit(1);
        }
    }
    else //cazul de inexistenta
    {
        printf("Programul nu va rula mai departe fara introducerea corecta a caii fisierului ce urmeaza a fi criptat!\n");
        exit(1);

    }

    //------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------//

    //partea a doua din enunt


    printf("\nIntrodu calea fisierului ce urmeaza sa fie folosita la pattern matching: ");
    char *image_name = read_fileBMP(); //citirea de la tastatura a caii imaginii pentru partea a 2-a a enuntului

    if(image_name != NULL) //caz existenta fisier
    {
        char image_name_grayscale[] = "test_grayscale.bmp";
        grayscale_image(image_name, image_name_grayscale);

        char *number_path;
        printf("\nCale pentru sablonul cifrei 0: ");
        number_path = read_fileBMP();
        if(number_path == NULL)
            exit(1); //oprirea algoritmului

        int i, j;
        unsigned int size_d = 0;

        window_pos *D = (window_pos *)malloc(sizeof(window_pos));
        window_pos *v;

        char *digit_name = (char *)malloc(120 * sizeof(char));

        //c este liniarizarea imaginii
        PIXEL *c = linearizeArray_BMP(image_name_grayscale);

        //colored este liniarizarea imaginii colorate
        PIXEL *colored = linearizeArray_BMP(image_name);

        //dimensiunea unui sablon citit de la tastatura
        FILE *f;
        f = fopen(number_path, "rb"); //nu exista cazul f == NULL deoarece stim sigur ca exista fisierul

        unsigned int number_widthW, number_heightH;
        fseek(f, 18, SEEK_SET);
        fread(&number_widthW, sizeof(unsigned int), 1, f);
        fread(&number_heightH, sizeof(unsigned int), 1, f);
        fclose(f);

        for(i = 0; i <= 9; i++)
        {
            unsigned int size = 0;
            PIXEL color;

            //schema de culori impusa impreuna cu citirea caii imaginilor
            if(i == 0)
                    {
                        strcpy(digit_name, number_path); //deja a fost citit, il pun direct fara citire de la tastatura
                        color.R = 255;
                        color.G = 0;
                        color.B = 0;
                    }
                else if(i == 1)
                        {

                            printf("\nCale pentru sablonul cifrei 1: ");
                            number_path = read_fileBMP();
                            if(number_path == NULL)
                                exit(1);
                            strcpy(digit_name, number_path);

                            color.R = 255;
                            color.G = 255;
                            color.B = 0;
                        }
                else if(i == 2)
                        {
                            printf("\nCale pentru sablonul cifrei 2: ");
                            number_path = read_fileBMP();
                            if(number_path == NULL)
                                exit(1);
                            strcpy(digit_name, number_path);

                            color.R = 0;
                            color.G = 255;
                            color.B = 0;
                        }
                else if(i == 3)
                        {

                            printf("\nCale pentru sablonul cifrei 3: ");
                            number_path = read_fileBMP();
                            if(number_path == NULL)
                                exit(1);
                            strcpy(digit_name, number_path);

                            color.R = 0;
                            color.G = 255;
                            color.B = 255;
                        }
                else if(i == 4)
                        {

                            printf("\nCale pentru sablonul cifrei 4: ");
                            number_path = read_fileBMP();
                            if(number_path == NULL)
                                exit(1);
                            strcpy(digit_name, number_path);

                            color.R = 255;
                            color.G = 0;
                            color.B = 255;
                        }
                else if(i == 5)
                        {

                            printf("\nCale pentru sablonul cifrei 5: ");
                            number_path = read_fileBMP();
                            if(number_path == NULL)
                                exit(1);
                            strcpy(digit_name, number_path);

                            color.R = 0;
                            color.G = 0;
                            color.B = 255;
                        }
                else if(i == 6)
                        {

                            printf("\nCale pentru sablonul cifrei 6: ");
                            number_path = read_fileBMP();
                            if(number_path == NULL)
                                exit(1);
                            strcpy(digit_name, number_path);

                            color.R = 192;
                            color.G = 192;
                            color.B = 192;
                        }
                else if(i == 7)
                        {
                            printf("\nCale pentru sablonul cifrei 7: ");
                            number_path = read_fileBMP();
                            if(number_path == NULL)
                                exit(1);
                            strcpy(digit_name, number_path);

                            color.R = 255;
                            color.G = 140;
                            color.B = 0;
                        }
                else if(i == 8)
                        {
                            printf("\nCale pentru sablonul cifrei 8: ");
                            number_path = read_fileBMP();
                            if(number_path == NULL)
                                exit(1);
                            strcpy(digit_name, number_path);

                            color.R = 128;
                            color.G = 0;
                            color.B = 128;
                        }
                else if(i == 9)
                        {

                            printf("\nCale pentru sablonul cifrei 9: ");
                            number_path = read_fileBMP();
                            if(number_path == NULL)
                                exit(1);
                            strcpy(digit_name, number_path);

                            color.R = 128;
                            color.G = 0;
                            color.B = 0;
                        }

            printf("Algoritmul lucreaza pentru cifra %d\n", i);
            digit_name[strlen(digit_name)] = '\0';

            //crearea unui fisier bmp pentru a aplica grayscale fara sa modifice sablonul initial
            char *cifra_noua = (char *)malloc(120 * sizeof(char));
            strcpy(cifra_noua, "cifra.bmp");

            //aplicarea functiei grayscale pentru sablonul nou creat
            grayscale_image(digit_name, cifra_noua);

            //cautarea sablonului cifra in imaginea grayscale
            v = template_match(c, image_name_grayscale, cifra_noua, 0.5, &size, &size_d, &D, color.R, color.G, color.B);
            free(cifra_noua);
        }

        //sortarea vectorului ce are toate detectiile celor 10 cifre
        qsort(D, size_d, sizeof(window_pos), cmp);

        //stergerea detectiilor ce au corelatia mica si au ferestre care se suprapun
        D = delete_allValues(D, number_widthW, number_heightH, &size_d);

        //colorarea fiecarei detectii ramase dupa stergerea celor necorespunzatoare
        colored = colorize_allWindows(colored, image_name_grayscale, digit_name, D, size_d);

        //scrierea in fisier a imaginii colorate
        writeFile_linearizedArray(image_name_grayscale, colored);

        //eliberarea memoriei si stergerea fisierelor utilizate de care nu mai avem nevoie
        free(c);
        free(v);
        free(D);
        free(number_path);
        free(digit_name);
        remove("cifra.bmp");
        remove("test_grayscale.bmp");
    }

    return 0;
}
