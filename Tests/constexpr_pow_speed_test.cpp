#include <opencv2/opencv.hpp>
#include <opencv2/core.hpp>
#include <opencv2/highgui.hpp>

#include "../sq_matrix.h"

template<typename T>
T microsec_cexprpow_stress_test(T x, T p, size_t repeats)
{
	auto t = new T(0);
	auto b = std::chrono::high_resolution_clock::now();
	for(size_t i = 0; i < repeats; ++i)
		*t += dixelu::utils::constexpr_pow<T>(x, p);
	auto e = std::chrono::high_resolution_clock::now();
	return T(std::chrono::duration_cast<std::chrono::microseconds>(e - b).count()) / repeats;
}

template<typename T>
T microsec_stdpow_stress_test(T x, T p, size_t repeats)
{
	auto t = new T(0);
	auto b = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < repeats; ++i)
		*t += std::pow<T, T>(x, p);
	auto e = std::chrono::high_resolution_clock::now();
	return T(std::chrono::duration_cast<std::chrono::microseconds>(e - b).count()) / repeats;
}

int main()
{
	int size = 750;
	double range = 2.54684;

	constexpr auto t = dixelu::utils::constexpr_pow<double>(0.156, 0.98);

	cv::Mat cspow(size, size, CV_64F);
	cv::Mat pow(size, size, CV_64F);
	cv::Mat diff(size, size, CV_64F);
	for(int i = 0; i < size; ++i)
	{
		for (int j = 0; j < size; ++j)
		{
			double x = (double(j) / size) * range;
			double y = (double(i) / size) * range;
			cspow.at<double>(i, j) = log10(microsec_cexprpow_stress_test<double>(x, y, 50)+1);
			pow.at<double>(i, j) = log10(microsec_stdpow_stress_test<double>(x, y, 50)+1);
		}
		std::cout << i << std::endl;
	}

	double cspow_max;
	double pow_max;
	double diff_min;
	double diff_max;

	diff = (cspow - pow) * (cspow - pow);

	cv::minMaxLoc(cspow, NULL, &cspow_max);
	cv::minMaxLoc(pow, NULL, &pow_max);
	cv::minMaxLoc(diff, &diff_min, &diff_max);

	cv::normalize(cspow, cspow, 0, 255, cv::NormTypes::NORM_MINMAX, CV_8U);
	cv::normalize(pow, pow, 0, 255, cv::NormTypes::NORM_MINMAX, CV_8U);
	cv::normalize(diff, diff, 0, 255, cv::NormTypes::NORM_MINMAX, CV_8U);

	cv::applyColorMap(cspow, cspow, cv::COLORMAP_JET);
	cv::applyColorMap(pow, pow, cv::COLORMAP_JET);
	cv::applyColorMap(diff, diff, cv::COLORMAP_JET);

	std::cout << "Max elements of: cspow: " << cspow_max << "; pow: " << pow_max << std::endl;
	std::cout << "Diff: [" << diff_min << "; " << diff_max << "]" << std::endl;
	cv::imshow("cspow", cspow);
	cv::imshow("pow", pow);
	cv::imshow("diff", diff);
	cv::waitKey(0);

	// conclusion: *it is on avarage too slow to be used anywhere outside of constexpr calculations*
}