---
layout: course
title: Biostats 285
---

# Biostats 285 - Advanced Bayesian Computing - Spring 2021

## Syllabus

### Course Topics

This 4-unit course covers recent developments in the Markov chain Monte Carlo (MCMC) literature with an emphasis on the challenges associated with both big data and high-dimensional statistical models.  The course starts with a quick review of the theoretical foundations of MCMC, Metropolis-Hastings and the Gibbs sampler.  Next, we cover Hamiltonian Monte Carlo and Riemannian Manifold Hamiltonian Monte Carlo.  The course then pivots to surrogate proposal methods and stochastic gradient MCMC and finishes by reviewing recent advances in non-reversible MCMC.

The first portion of this course includes lectures and textbook reading assignments. After this, assignments include the reading of scholarly papers, and the majority of class time will be spent discussing these papers.  Although practical implementation is not emphasized, a final project involves the rigorous implementation of one or more of the methods encountered during the quarter.


### Lecture and Practicals

Tues, 1-230pm @ Zoom <https://ucla.zoom.us/j/7120566052>\
Thur, 1-230pm @ Zoom <https://ucla.zoom.us/j/7120566052>

### Instructor

[Andrew J. Holbrook](http://andrewjholbrook.github.io/)\
Email: <aholbroo@g.ucla.edu>\
Office hours: offered generously; please email me.

### Course Webpage

<https://ucla-biostats-285.github.io>

### Course Learning Objectives and Competencies

1. Know and understand core MCMC algorithms
2. Understand Metropolis-Hastings-Green algorithm
3. Gain familiarity with various adaptive MCMC strategies
4. Know and understand Hamiltonian Monte Carlo and extensions
5. Understand MCMC challenges associated with big data and large models
6. Gain familiarity with various scalable MCMC approaches
7. Know and understand basic theory of nonreversible MCMC
8. Understand tradeoffs between different MCMC strategies
9. Understand how different MCMC strategies may be combined

### Prerequisite

Students are expected to have a firm understanding of Bayesian statistics and the basics of Markov chain Monte Carlo, having passed courses on the levels of Biostats M234 and 202C.  Moreover, successful completion of the final project will require fluency in statistical programming (in, e.g., R, python, Julia, etc.), having passed a course on the levels of Biostats 203A and 213.

### Course Texts

Since the course will cover a collection of specialized topics in a variety of areas, there will not be a single text. Instead, chapters and topics will be sampled from multiple different texts. The topics covered will roughly follow material at the level of the following books and papers:

1. Ronald Christensen, Wesley O. Johnson, Adam J. Branscum, Timothy E. Hanson (2010) Bayesian Ideas and Data Analysis: An Introduction for Scientists and Statisticians: Chapter 6

2. Brooks, Steve, Andrew Gelman, Galin Jones, and Xiao-Li Meng, eds. Handbook of Markov Chain Monte Carlo. CRC press, 2011: Chapters 1, 3, 4 and 5.

3. Girolami, Mark, and Ben Calderhead. "Riemann manifold langevin and hamiltonian monte carlo methods." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 73.2 (2011): 123-214.

4. Zhang, Cheng, Babak Shahbaba, and Hongkai Zhao. "Hamiltonian Monte Carlo acceleration using surrogate functions with random bases." Statistics and computing 27, no. 6 (2017): 1473-1490.

5. Levy, Daniel, Matt D. Hoffman, and Jascha Sohl-Dickstein. "Generalizing Hamiltonian Monte Carlo with Neural Networks." International Conference on Learning Representations. 2018.

6. Nemeth, Christopher, and Paul Fearnhead. "Stochastic gradient markov chain monte carlo." Journal of the American Statistical Association (2020): 1-47.

7. Bouchard-Côté, Alexandre, Sebastian J. Vollmer, and Arnaud Doucet. "The bouncy particle sampler: A nonreversible rejection-free Markov chain Monte Carlo method." Journal of the American Statistical Association 113.522 (2018): 855-867.

8. Fearnhead, Paul, et al. "Piecewise deterministic Markov processes for continuous-time Monte Carlo." Statistical Science 33.3 (2018): 386-412.

9. Zanella, Giacomo. "Informed proposals for local MCMC in discrete spaces." Journal of the American Statistical Association 115.530 (2020): 852-865.

10. Gagnon, Philippe, and Arnaud Doucet. "Nonreversible Jump Algorithms for Bayesian Nested Model Selection." Journal of Computational and Graphical Statistics (2020): 1-12.




### Grading

Grades will be based on the homework (50%) and final phylogenetic analysis project (50%).

### Project

A Final Project consists of either (1) the replication of results and comparison between two or more methods encountered in the course, (2) the development and successful application of a prototype method inspired by course content or (3) the replication of results and comparison between two or more methods to be agreed upon between student and instructor. One page project proposals are due 5/13 by 5PM via email (PDF). Statistical code and project manuscript (max. 10 typed pages) must be submitted to me via email (PDF) by the final day of instruction (6/3) by 5PM.

### Attendance

I usually expect class attendance. Given the COVID-19 pandemic, I understand students may not be able to attend certain lectures due to medical or technical reasons. Please proactively communicate with me about your circumstances at earliest chance.

### UCLA ADA Policy

Students needing academic accommodation based on a disability should contact the Center for Accessible Education (CAE) at (310)825-1501 or in person at Murphy Hall A255. When possible, students should contact the CAE within the first two weeks of the term as reasonable notice is needed to coordinate accommodations. For more information visit <https://www.cae.ucla.edu>.

ADA Contact:
Nickey Woods
Center for Accessible Education
A255 Murphy Hall.
Phone: (310)825-1501
TTY/TTD: (310)206-6083
Fax: (310)825-9656

### Inclusivity

UCLA’s Office for Equity, Diversity, and Inclusion provides resources, events, and information about current initiatives at UCLA to support equality for all members of the UCLA community. I hope that you will communicate with me or your TA if you experience anything in this course that does not support an inclusive environment, and you can also report any incidents you may witness or experience on campus to the Office of Equity, Diversity, and Inclusion on their website <https://equity.ucla.edu>.


### Academic Integrity

**Message about Academic Integrity to All UCLA Students from UCLA Dean of Students**: UCLA is a community of scholars. In this community, all members including faculty staff and students alike are responsible for maintaining standards of academic honesty. As a student and member of the University community, you are here to get an education and are, therefore, expected to demonstrate integrity in your academic endeavors. You are evaluated on your own merits. Cheating, plagiarism, collaborative work, multiple submissions without the permission of the professor, or other kinds of academic dishonesty are considered unacceptable behavior and will result in formal disciplinary proceedings usually resulting in suspension or dismissal.

**Forms of Academic Dishonesty**: As specified in the UCLA Student Conduct Code, violations or attempted violations of academic dishonesty include, but are not limited to, cheating, fabrication, plagiarism, multiple submissions or facilitating academic integrity:
• Allowing another person to take a quiz, exam, or similar evalution for you
• Using unauthorized material, information, or study aids in any academic exercise or examination – textbook, notes, formula list, calculators, etc.
• Unauthorized collaboration in providing or requesting assistance, such as sharing information
• Unauthorized use of someone else’s data in completing a computer exercise
• Altering a graded exam or assignment and requesting that I be regraded

**Plagiarism**: Presenting another’s words or ideas as if they were one’s own
• Submitting as your own through purchase or otherwise, part of or an entire work produced verbatim by someone else
• Paraphrasing ideas, data or writing without properly acknowledging the source
• Unauthorized transfer and use of someone else’s computer file as your own
• Unauthorized use of someone else’s data in completing a computer exercise

**Multiple Submissions**: Submitting the same work (with exact or similar content) in more than one class without permission from the instructor to do so. This includes courses you are currently taking, as well as courses you might take in another quarter.

**Facilitating Academic Dishonesty**: Participating in any action that compromises the integrity of the academic standards of the University; assisting another to commit an act of academic dishonesty
• Taking a quiz, exam, or similar evaluation in place of another person
• Allowing another student to copy from you
• Providing material or other information to another student with knowledge that such assistance could be used in any of the violations stated above (e.g., giving test information to students in other discussion sections of the same course)
• Altering data to support research
• Presenting results from research that was not performed
• Crediting source material that was not used for  research

While you are here at UCLA, if you are unsure whether what you are considering doing is cheating, **don’t take chances** – ask your professor. In addition, avoid placing yourself in situations which might lead your professor to **suspect you of cheating**.

**Alternatives to Academic Dishonesty**

• Seek out help – Meet with your professor, ask for assistance as needed.
• Ask for an extension – if you explain your situation to your professor, she/he might be able to grant you an extended deadline for an upcoming assignment.
•	See a counselor at Student Psychological Services, and/or your school, college or department – UCLA has many resources for students who are feeling the stresses of academic and personal pressures.

If you would like more information, please come see us at the Dean of Students’ Office in 1206 Murphy Hall, call us at (310)825-3871 or visit their website at <https://www.deanofstudents.ucla.edu>.
