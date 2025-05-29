#ifndef _EXPORT_
#define _EXPORT_

#define TINYEXR_USE_MINIZ 0
#define TINYEXR_USE_STB_ZLIB 1 
#define TINYEXR_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define STB_IMAGE_WRITE_STATIC
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image_write.h"
#include "stb_image.h"
#include "tinyexr.h"

// export data to C
void writeTabC(mat3 * tab, vec2 * tabMagFresnel, int N)
{
    ofstream file("results/ltc.inc");

    file << std::fixed;
    file << std::setprecision(6);

    file << "static const int size = " << N  << ";" << endl << endl;

    file << "static const mat33 tabM[size*size] = {" << endl;
    for (int t = 0; t < N; ++t)
    for (int a = 0; a < N; ++a)
    {
        file << "{";
        file << tab[a + t*N][0][0] << ", " << tab[a + t*N][0][1] << ", " << tab[a + t*N][0][2] << ", ";
        file << tab[a + t*N][1][0] << ", " << tab[a + t*N][1][1] << ", " << tab[a + t*N][1][2] << ", ";
        file << tab[a + t*N][2][0] << ", " << tab[a + t*N][2][1] << ", " << tab[a + t*N][2][2] << "}";
        if (a != N - 1 || t != N - 1)
            file << ", ";
        file << endl;
    }
    file << "};" << endl << endl;

    file << "static const mat33 tabMinv[size*size] = {" << endl;
    for (int t = 0; t < N; ++t)
    for (int a = 0; a < N; ++a)
    {
        mat3 Minv = glm::inverse(tab[a + t*N]);

        file << "{";
        file << Minv[0][0] << ", " << Minv[0][1] << ", " << Minv[0][2] << ", ";
        file << Minv[1][0] << ", " << Minv[1][1] << ", " << Minv[1][2] << ", ";
        file << Minv[2][0] << ", " << Minv[2][1] << ", " << Minv[2][2] << "}";
        if (a != N - 1 || t != N - 1)
            file << ", ";
        file << endl;
    }
    file << "};" << endl << endl;

    file << "static const float tabMagnitude[size*size] = {" << endl;
    for (int t = 0; t < N; ++t)
    for (int a = 0; a < N; ++a)
    {
        file << tabMagFresnel[a + t*N][0] << "f";
        if (a != N - 1 || t != N - 1)
            file << ", ";
        file << endl;
    }
    file << "};" << endl;

    file.close();
}

// export data to MATLAB
void writeTabMatlab(mat3 * tab, vec2 * tabMagFresnel, int N)
{
    ofstream file("results/ltc.mat");

    file << "# name: tabMagnitude" << endl;
    file << "# type: matrix" << endl;
    file << "# ndims: 2" << endl;
    file << " " << N << " " << N << endl;

    for (int t = 0; t < N; ++t)
    {
        for (int a = 0; a < N; ++a)
            file << tabMagFresnel[a + t*N][0] << " ";
        file << endl;
    }

    for (int row = 0; row < 3; ++row)
    for (int column = 0; column < 3; ++column)
    {

        file << "# name: tab" << column << row << endl;
        file << "# type: matrix" << endl;
        file << "# ndims: 2" << endl;
        file << " " << N << " " << N << endl;

        for (int t = 0; t < N; ++t)
        {
            for (int a = 0; a < N; ++a)
                file << tab[a + t*N][column][row] << " ";
            file << endl;
        }

        file << endl;
    }

    file.close();
}

// export data to DDS
#include "dds.h"
#include "float_to_half.h"

void writeDDS(const char* path, float* data, int N)
{
    EXRImage image;
    InitEXRImage(&image);

    image.num_channels = 4;

    int width = N, height = N;

    std::vector<float> images[4];
    images[0].resize(width * height); // R
    images[1].resize(width * height); // G
    images[2].resize(width * height); // B
    images[3].resize(width * height); // A

    // OpenEXR 要求每个通道为独立数组，按 R、G、B、A 分离
    for (int i = 0; i < width * height; ++i) {
        images[0][i] = data[4 * i + 0];
        images[1][i] = data[4 * i + 1];
        images[2][i] = data[4 * i + 2];
        images[3][i] = data[4 * i + 3];
    }

    float* image_ptr[4];
    image_ptr[0] = images[0].data();
    image_ptr[1] = images[1].data();
    image_ptr[2] = images[2].data();
    image_ptr[3] = images[3].data();

    image.images = (unsigned char**)image_ptr;
    image.width = width;
    image.height = height;

    EXRHeader header;
    InitEXRHeader(&header);

    header.num_channels = 4;
    header.channels = (EXRChannelInfo*)malloc(sizeof(EXRChannelInfo) * 4);
    strncpy(header.channels[0].name, "R", 255);
    strncpy(header.channels[1].name, "G", 255);
    strncpy(header.channels[2].name, "B", 255);
    strncpy(header.channels[3].name, "A", 255);

    // 所有通道使用 HALF (16-bit float)
    header.pixel_types = (int*)malloc(sizeof(int) * 4);
    header.requested_pixel_types = (int*)malloc(sizeof(int) * 4);
    for (int i = 0; i < 4; i++) {
        header.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT; // 输入数据类型
        header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_HALF; // 输出保存格式
    }

    const char* err = nullptr;
    int ret = SaveEXRImageToFile(&image, &header, path, &err);
    if (ret != TINYEXR_SUCCESS) {
        std::cerr << "Failed to save EXR: " << err << std::endl;
        FreeEXRErrorMessage(err);
    } else {
        std::cout << "Saved EXR" << std::endl;
    }

    free(header.channels);
    free(header.pixel_types);
    free(header.requested_pixel_types);
}

void writeDDS(vec4* data1, vec4* data2, int N)
{
    writeDDS("results/ltc_1.exr", &data1[0][0], N);
    writeDDS("results/ltc_2.exr", &data2[0][0], N);
}

// export data to Javascript
void writeJS(vec4* data1, vec4* data2, int N)
{
    ofstream file("results/ltc.js");

    file << "var g_ltc_1 = [" << endl;

    for (int i = 0; i < N*N; ++i)
    {
        // store the variable terms
        file << data1[i].x << ", ";
        file << data1[i].y << ", ";
        file << data1[i].z << ", ";
        file << data1[i].w << ", " << endl;
    }
    file << "];" << endl;

    file << "var g_ltc_2 = [";
    for (int i = 0; i < N*N; ++i)
    {
        file << data2[i].x << ", ";
        file << data2[i].y << ", ";
        file << data2[i].z << ", ";
        file << data2[i].w << ", " << endl;
    }
    file << "];" << endl;

    file.close();
}

#endif
