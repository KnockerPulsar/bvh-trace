#include "vec_math.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

class Image {
public:
  std::vector<char> pixels;
  int               width, height;

  Image() = default;

  Image(std::string path) {
    std::ifstream in(path, std::ios::binary);

    std::string s;
    in >> s;

    if (s != "P6") {
      exit(1);
    }

    // Skip comments at the start of the ppm file
    for (;;) {
      std::getline(in, s);
      if (s.empty()) {
        continue;
      }

      if (s[0] != '#') {
        break;
      }
    }

    std::stringstream str(s);
    int               width;
    int               height;
    int               maxColor;

    str >> width >> height;
    in >> maxColor;

    if (maxColor != 255) {
      exit(1);
    }

    // Skip until the end of the line?
    {
      std::string tmp;
      std::getline(in, tmp);
    }

    std::vector<char> data(width * height * 3);
    in.read(reinterpret_cast<char *>(data.data()), data.size());

    pixels = data;
    this->width  = width;
    this->height = height;
  }

  static Image magenta(int width, int height) {
    Image result;
    result.width = width;
    result.height = height;

    for (int i = 0; i < width*height; i++) {
      result.pixels.push_back((char)255);
      result.pixels.push_back((char)0);
      result.pixels.push_back((char)255);
    }

    return result;
  }

  static Image fromColor(int width, int height, float3 color) {
    Image result;
    result.width = width;
    result.height = height;

    for (int i = 0; i < width*height; i++) {
      result.pixels.push_back((char)color.x);
      result.pixels.push_back((char)color.y);
      result.pixels.push_back((char)color.z);
    }

    return result;
  }

  void save(std::string path) {
    std::ofstream out(path, std::ios::binary);

    out << "P6\n";
    out << width << " " << height << "\n";
    out << "255\n";

    out.write(pixels.data(), pixels.size());
  }

  Image toRGBA() {
    Image result;
    result.width  = width;
    result.height = height;

    for (int i = 0; i < pixels.size(); i += 3) {
      result.pixels.push_back(pixels[i + 0]);
      result.pixels.push_back(pixels[i + 1]);
      result.pixels.push_back(pixels[i + 2]);
      result.pixels.push_back(0);
    }

    return result;
  }

  Image toRGB() {
    Image result;
    result.width  = width;
    result.height = height;

    for (int i = 0; i < pixels.size(); i += 4) {
      result.pixels.push_back(pixels[i + 0]);
      result.pixels.push_back(pixels[i + 1]);
      result.pixels.push_back(pixels[i + 2]);
    }

    return result;
  }

  void writeToPixel(int x, int y, char value) {
    int pixelIndex = (x + y * width) * 3;
    pixels[pixelIndex + 0] = value;
    pixels[pixelIndex + 1] = value;
    pixels[pixelIndex + 2] = value;
  }

  void writeToPixel(int x, int y, float3 color) {
    int pixelIndex = (x + y * width) * 3;
    pixels[pixelIndex + 0] = color.x;
    pixels[pixelIndex + 1] = color.y;
    pixels[pixelIndex + 2] = color.z;
  }
};
