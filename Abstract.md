A Greedy Planning Method for Multi-Drone Fleets Control

Currently multi-drone fleets operating a variety of missions in common airspace are undergoing rapid development, which will require trajectory planning to fulfill the mission requirements as well as guarantee the safety of the drone systems. Recognizing this, a recently published paper [1] provided a Fly-by-Logic method to solve this problem of safe planning and control for multi-agent systems across multiple tasks. In accordance with the Fly-by-Logic, Signal Temporal Logic (STL) is being used for mission specification language and generating waypoints by maximizing the robustness the STL mission.

	It is necessary for the Fly-by-Logic method to specify all the mission requirements before planning, which is a centralized method that plan for the whole system all at once. This method could not effectively react to emergencies or random variations under the complicated airspace circumstances due to the situation that a slight change of one single drone in the fleet would cause the system to be planned all over again. A sequential planning method which plans in serial based on the Fly-by-Logic method would enable each drone to alter its trajectory with minimized influence imposing on other drones. Business concerns on security and privacy could be eliminated with the information exchange narrowing down to single drone but not the entire drone fleet.

	Sequential planning will separate the mission specifications and rank the missions with pre-set priority. The mission of the first drone in sequence will be planned first and have the maximum free airspace. The next drone will be planned with trajectories of the previous drone considered as obstacles in the airspace. Due to the limited airspace the latter drones can be planned in, the robustness of the multi-drone system under sequential planning will be slightly less than the one under centralized planning. However, the robustness is positive, which means this multi-agent system is safe, for almost 95% cases and the exceptional cases happen when the airspace is extremely congested. Considering the situation that most airspace is currently unexploited, this sequential planning method is suitable for current application.

	The planning time of the sequential planning decrease dramatically compared to the centralized one. Even though the centralized method only plans once for the whole system, the computing complexity is considerably greater than the sequential one, which reflects to the planning time. With the number of drones increasing in the drone fleet, the planning time of sequential method can reach up to 4 times faster than the centralized one. The sequential method will provide a fast response to external disturbance and efficient information interaction thus guarantee safety to some extent. The large-scale drone fleet typically with more uncertainty can benefit more from this sequential method.