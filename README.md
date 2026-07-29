<div style="max-width: 600px; margin: auto;">

# irspec - finishing up and writing down some end goals

With the help of chat we have made more progress than otherwise, but we have several threads up in the air that we need to make sure we at least mark before we move on to school starting and the preparations around that. 

## hitherto

These were the problems leading into this summer that we knew about and wanted to work on.

1. Joe wrote some code that was slightly wrong for difficult to understand reasons.
2. Milan's code was difficult to understand what it was doing since it was someone else's code and also written in Matlab.
3. Milan also produced a set of index of refraction and absorption (extinction) coefficients and we had no idea how he did this. We were told it was based off a paper (Herguedas) and that there was some polynomial fitting involved to get the index of refraction at other concentrations of O:Si. 
4. Kramers-Kronig relations are something that I genuinely want to understand and so we need to pull it off and then play with it.
5. Rocha's code is not completely correct and also difficult to understand. We need to work to improve it and adapt it to our case of reflection.

## herein

This is what we worked on and where it is at this moment as we are coming to the end of the summer

1. Joe's code we corrected. This took some time since we had to get into what exactly he was doing and what the issues were. We also had Allison's experimental results and computational results. By comparing these we realized that the issue was focused on our definition of Absorption. It turns out that for reflection we should have been referencing the case of gold only with no film on top. BUT WE STILL NEED TO VERIFY THIS. This seems like the right idea and an improvement to our code, and it was what Milan was doing as well, but that took a good bit of digging to figure out.
2. Milan's code we had chat translate to python. The crucial bit here was verifying how exactly he was calculating the Absorption and then we were able to proceed with our matrix method and reproduce similar plots.
3. Crucial to the above step was coming up with a similar method for producing the index of refraction at any concentration of O:Si that we wanted. Some of the concentrations that had been identified as being between the ones that we had n and k generated for, so we needed to be able to reproduce this. This also gave us a better understanding of the data that this came from and an understanding of what error's it may introduce.
4. Kramers-Kronig is much harder than I anticipated, but I should have anticipated it being tough, because it is an integral over an infinite range. That means that experimentally or computationally there are always some models or built in assumptions that are at play. But we worked on some methods for doing this and that is ongoing at this moment.
5. The Rocha code was finally found and we are currently working on adapting it to our needs. One thing we learned is that it is only for transmission, not reflection, so we need to add that capability. For another, it only accepts a single value for the index of refraction of the substrate, so we need to make it capable of doing this over a range where the index changes by significant amounts as it does in gold. The goal with this piece of software is to extract the n and k values from an experimental absorption spectra where we do not have an underlying model of what the n and k values should be.

## hence

1. This code now works, but we need to verify with some citations these code changes. Then, other tasks follow below.
2. Now that we have a translation of this, we need to see where exactly this method differs from our own matrix method, and decide on why and which one is better.
3. Plot the absorption as a function of concentration and verify that that changes the peak position, and that changing the thickness `d` changes the height. This is a claim made by Gregory that I would like to see verified. Make sure we have a function where any concentration of O:Si can be used and the data can be saved or plotted.
4. I want to be able to put *any* function into this and take the Kramers-Kronig of it. Whether it is realistic or not. So I want a library of different functions over which the Kramers-Kronig can be taken, and I want the ability to make my own function and then KK it.
5. We need to adapt this to reflection and make it a complex list of index of refraction and absorption coefficients. We currently have some corrections already in place, but they are very chatgpt heavy. We need to proceed with caution but proceed nonetheless and work on validating these functions and expanding the role of this software. This currently needs to be expanded in several ways:
   a. the equations need to be written for reflection from a substrate
   b. the equations need to have a changeable angle of incidence since it assumes normal incidence at the moment
   c. the equations need to be corrected since there are some conceptual issues with the original script
   d.



    [NbConvertApp] Converting notebook README.ipynb to markdown
    [NbConvertApp] Writing 5222 bytes to README.md


</div>
