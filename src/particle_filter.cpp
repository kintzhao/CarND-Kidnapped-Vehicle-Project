/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    num_particles = 100;
    std::normal_distribution<double> dist_x(x, std[0]);
    std::normal_distribution<double> dist_y(y, std[1]);
    std::normal_distribution<double> dist_theta(theta, std[2]);
    particles.resize(num_particles);
    std::cout<<"init"<<std::endl;
    for (size_t i = 0; i < num_particles; i++) {
        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
        particles[i].id = i;
        particles[i].weight = 1;
        //std::cout<<i<<": x y theta"<<particles[i].x <<" "<<particles[i].y <<" "<<particles[i].theta<<std::endl;
    }
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    std::normal_distribution<double> N_x(0, std_pos[0]);
    std::normal_distribution<double> N_y(0, std_pos[1]);
    std::normal_distribution<double> N_theta(0, std_pos[2]);
    std::cout<<"prediction"<<std::endl;
    for(size_t i = 0; i < num_particles; i++){
        //std::cout<<i<<": x y theta"<<particles[i].x <<" "<<particles[i].y <<" "<<particles[i].theta<<std::endl;
        if( fabs(yaw_rate) < 0.0001){
            particles[i].x += velocity * delta_t * cos(particles[i].theta);
            particles[i].y += velocity * delta_t * sin(particles[i].theta);

        } else{
            particles[i].x += velocity / yaw_rate * ( sin( particles[i].theta + yaw_rate * delta_t ) - sin(particles[i].theta) );
            particles[i].y += velocity / yaw_rate * ( cos( particles[i].theta ) - cos( particles[i].theta + yaw_rate * delta_t) );
            particles[i].theta += yaw_rate * delta_t;
        }

        particles[i].x += N_x(gen);
        particles[i].y += N_y(gen);
        particles[i].theta += N_theta(gen);
        // particles[i].id = i;
        //std::cout<<"=>predict: "<<i<<": x y theta"<<particles[i].x <<" "<<particles[i].y <<" "<<particles[i].theta<<std::endl;
    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    for (auto& obs: observations) {
        double minDist = std::numeric_limits<double>::max();
        for (auto pred: predicted) {
            double distance = dist(pred.x, pred.y, obs.x, obs.y);
            if (distance < minDist) {
                minDist = distance;
                obs.id = pred.id;
            }
            //std::cout<<" minDist:"<<std::numeric_limits<double>::max()<<std::endl;
        }
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    for (auto& p: particles) {
        vector<LandmarkObs> predictions;
        for (auto landmark: map_landmarks.landmark_list) {
            if (dist(p.x, p.y, landmark.x_f, landmark.y_f) < sensor_range) {
                predictions.push_back(LandmarkObs{landmark.id_i, landmark.x_f, landmark.y_f});
            }
        }

        vector<LandmarkObs> observation_transformed;
        for(auto obs: observations){
            LandmarkObs tmp;
            tmp.x = obs.x * cos(p.theta) - obs.y * sin(p.theta) + p.x;
            tmp.y = obs.x * sin(p.theta) + obs.y * cos(p.theta) + p.y;
            tmp.id = obs.id;
            observation_transformed.push_back(tmp);
        }

        dataAssociation(predictions, observation_transformed);

        p.weight = 1.0;
        for(auto obs_t: observation_transformed){
            Map::single_landmark_s landmark = map_landmarks.landmark_list.at(obs_t.id - 1);
            double x_term = pow(obs_t.x - landmark.x_f, 2) / (2 * pow(std_landmark[0], 2));
            double y_term = pow(obs_t.y - landmark.y_f, 2) / (2 * pow(std_landmark[1], 2));
            double w = exp(-(x_term + y_term)) / (2 * M_PI * std_landmark[0] * std_landmark[1]);
            p.weight *=  w;
            cout<<"x_term = "<<x_term<<endl;
            cout<<"y_term = "<<y_term<<endl;
            cout<<"w = "<<w<<endl;
        }

        weights.push_back(p.weight);
    }
//    for (int i = 0; i < weights.size(); i++) {
//        cout<<"weight[" <<i<<"] = "<<weights[i]<<endl;
//    }


}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    vector<Particle> resampled_particles;
    double max_weight = *max_element(weights.begin(), weights.end());
    std::normal_distribution<double> dist(0, 2 * max_weight);
    int idx = 0;
    for (size_t i = 0; i < num_particles; i++) {
        double beta = dist(gen);
        while (beta > weights[idx]) {
            beta -= weights[idx];
            idx = (idx + 1) % num_particles;
        }
        resampled_particles.push_back(particles[idx]);
    }
    particles = resampled_particles;
    weights.clear();
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
