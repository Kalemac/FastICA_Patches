module FastICA_Patches {
	exports org.fastica.swing;
	exports org.fastica.math;
	exports org.fastica;
	exports org.fastica.util;
	exports FastICA_Patches;

	requires commons.math3;
	requires ejml.core;
	requires java.desktop;
	requires java.prefs;
	requires opencsv;
}